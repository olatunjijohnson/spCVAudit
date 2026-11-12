.cov_matrix <- function(x, y, model, same = FALSE) {
  d <- .pair_distance(x, y)
  k <- matrix(model$variance * spatial_correlation(model, as.vector(d)), nrow(x), nrow(y))
  if (same) diag(k) <- model$variance
  k
}

.gaussian_reduction <- function(x, q, model) {
  s <- .cov_matrix(x, x, model, same = TRUE)
  b <- .cov_matrix(q, x, model)
  u <- tryCatch(chol(s), error = function(e)
    stop("Training covariance is not positive definite; check coincident sites and the nugget.", call. = FALSE))
  v <- forwardsolve(t(u), t(b))
  explained <- colSums(v^2)
  distance <- .pair_distance(q, x)
  first <- max.col(-distance, ties.method = "first")
  b1 <- b[cbind(seq_len(nrow(q)), first)]
  nearest <- b1^2 / model$variance
  additional <- explained - nearest
  tol <- 1e-9 * model$variance
  if (any(additional < -tol) || any(explained > model$variance + tol))
    stop("Numerical covariance calculations are unstable for these parameters.", call. = FALSE)
  additional <- pmax(additional, 0)
  if (model$nugget > 0) {
    residual <- b - t(s[, first, drop = FALSE]) * (b1 / model$variance)
    bound <- pmax(0, pmin(rowSums(residual^2) / model$nugget,
                          model$psill - nearest))
  } else bound <- rep(NA_real_, nrow(q))
  list(mse = pmax(model$variance - explained, 0), explained = explained,
       nearest = nearest, additional = additional, bound = bound,
       nearest_distance = distance[cbind(seq_len(nrow(q)), first)],
       weights = backsolve(u, v), conditional_variance = pmax(model$variance - explained, 0))
}

#' Calculate Gaussian reference MSE, its decomposition, and a bound
#'
#' Assumes jointly Gaussian responses with known mean and the supplied stationary
#' covariance, fixed locations, folds and target weights. Prediction uses all
#' available observations. Target observations have fresh independent noise even
#' when their coordinates coincide with training locations. The calculation
#' averages over repeated response datasets; it is not an observed CV score.
#'
#' `optimism` is target MSE minus expected CV MSE. It equals `one_neighbour`
#' (total variance times `J2`) plus `additional_neighbours`. The latter is the
#' validation-minus-target difference in the extra MSE reduction supplied by
#' observations beyond the nearest. With strictly positive nugget, `lower` and
#' `upper` bound optimism using conditional covariance residuals. These are
#' deterministic limits under the supplied model, not confidence intervals.
#' Bounds are `NA` when the nugget is zero. Supplying estimated covariance gives
#' a plug-in model calculation, not a guarantee about actual CV bias.
#'
#' @param locations Numeric two-column coordinates or projected spatial points.
#' @param folds A partition accepted by [cv_folds()].
#' @param target A [prediction_target()] or target coordinates.
#' @param covariance A [spatial_covariance()] object. A correlation function
#'   alone is insufficient because the calculation requires total variance.
#' @return A `spcv_gaussian` list with scalar `target_mse`, `cv_mse`, `optimism`,
#'   `J2`, `one_neighbour`, `additional_neighbours`, `lower`, `upper`, and
#'   `estimated_covariance`; and location-level `validation` and `target` data
#'   frames of MSE reductions and bounds. Costs are cubic in training-set size.
#' @export
#' @examples
#' xy <- cbind(c(-1, 0, 1), 0)
#' q <- prediction_target(matrix(c(2, 0), nrow = 1))
#' result <- gaussian_cv(xy, 1:3, q, spatial_covariance(1, 0.8))
#' result
#' # Identical nearest distances, but full-neighbour optimism is about 0.02306.
gaussian_cv <- function(locations, folds, target, covariance) {
  if (!inherits(covariance, "spcv_covariance"))
    stop("`covariance` must be a spatial_covariance() object.", call. = FALSE)
  z <- .target_coordinates(locations, target)
  n <- nrow(z$x)
  folds <- cv_folds(folds, n)
  columns <- c("mse", "explained", "nearest", "additional", "bound", "nearest_distance")
  cv <- as.data.frame(matrix(NA_real_, n, length(columns), dimnames = list(NULL, columns)))
  for (f in folds) {
    r <- .gaussian_reduction(z$x[f$train, , drop = FALSE], z$x[f$test, , drop = FALSE], covariance)
    for (nm in columns) cv[f$test, nm] <- r[[nm]]
  }
  r <- .gaussian_reduction(z$x, z$q, covariance)
  qt <- as.data.frame(r[columns])
  w <- z$target$weights
  j2 <- mean(cv$nearest) / covariance$variance - sum(w * qt$nearest) / covariance$variance
  one <- covariance$variance * j2
  additional <- mean(cv$additional) - sum(w * qt$additional)
  target_mse <- sum(w * qt$mse)
  cv_mse <- mean(cv$mse)
  structure(list(target_mse = target_mse, cv_mse = cv_mse, optimism = target_mse - cv_mse,
                 J2 = j2, one_neighbour = one, additional_neighbours = additional,
                 lower = one - sum(w * qt$bound), upper = one + mean(cv$bound),
                 estimated_covariance = isTRUE(covariance$estimated),
                 validation = cv, target = qt, target_weights = w), class = "spcv_gaussian")
}

#' @export
print.spcv_gaussian <- function(x, ...) {
  cat("Gaussian reference calculation", if (x$estimated_covariance) " (estimated covariance; plug-in)", "\n", sep = "")
  print(unlist(x[c("target_mse", "cv_mse", "optimism", "one_neighbour",
                   "additional_neighbours", "lower", "upper")]))
  invisible(x)
}

#' Evaluate conditional target MSE under a specified Gaussian distribution
#'
#' This reference calculation uses the supplied mean and covariance as the
#' generating distribution. Conditional on observed training responses, target
#' MSE is conditional variance plus squared difference between the supplied
#' predictions and the true conditional mean. It is available in a simulation
#' with known generating parameters; estimated parameters do not make it an
#' exact evaluation of real-world prediction error.
#'
#' @param locations Numeric two-column training coordinates or projected points.
#' @param response Finite training responses.
#' @param target A [prediction_target()] or target coordinates.
#' @param covariance The specified generating [spatial_covariance()].
#' @param predictions Finite predictions at each target location. If `NULL`,
#'   uses the conditional mean under the supplied Gaussian distribution.
#' @param mean Known constant mean of the generating distribution.
#' @return A list with weighted `mse`, vectors `pointwise_mse`, `conditional_mean`,
#'   `conditional_variance`, and `predictions`.
#' @export
#' @examples
#' xy <- cbind(c(0, 1, 2), 0)
#' gaussian_target_mse(xy, c(0.2, -0.1, 0.5), matrix(c(1.5, 0), 1),
#'                     spatial_covariance(1, 0.8))
gaussian_target_mse <- function(locations, response, target, covariance,
                                 predictions = NULL, mean = 0) {
  if (!inherits(covariance, "spcv_covariance"))
    stop("`covariance` must be a spatial_covariance() object.", call. = FALSE)
  z <- .target_coordinates(locations, target)
  .scalar(mean, "mean")
  if (!is.numeric(response) || length(response) != nrow(z$x) || any(!is.finite(response)))
    stop("`response` must have one finite value per training observation.", call. = FALSE)
  r <- .gaussian_reduction(z$x, z$q, covariance)
  cm <- as.numeric(mean + crossprod(r$weights, response - mean))
  if (is.null(predictions)) predictions <- cm
  if (!is.numeric(predictions) || length(predictions) != nrow(z$q) || any(!is.finite(predictions)))
    stop("`predictions` must have one finite value per target location.", call. = FALSE)
  point <- r$conditional_variance + (cm - predictions)^2
  list(mse = sum(z$target$weights * point), pointwise_mse = point,
       conditional_mean = cm, conditional_variance = r$conditional_variance,
       predictions = predictions)
}
