.cdf_areas <- function(x, y, wx, wy) {
  cuts <- sort(unique(c(x, y)))
  if (length(cuts) < 2L) return(c(positive = 0, negative = 0, W = 0))
  mass <- numeric(length(cuts))
  for (i in seq_along(x)) mass[match(x[i], cuts)] <- mass[match(x[i], cuts)] + wx[i]
  for (i in seq_along(y)) mass[match(y[i], cuts)] <- mass[match(y[i], cuts)] - wy[i]
  difference <- utils::head(cumsum(mass), -1L)
  widths <- diff(cuts)
  p <- sum(widths * pmax(difference, 0))
  n <- sum(widths * pmax(-difference, 0))
  c(positive = p, negative = n, W = p + n)
}

#' Compare validation and target nearest-neighbour distances
#'
#' Calculates exact integrals for weighted empirical distributions, without a
#' numerical integration grid. `A` is mean target distance minus mean validation
#' distance; `W` is their one-dimensional Wasserstein distance. `J1` and `J2` are
#' validation minus target averages of correlation and squared correlation.
#' `U1` and `U2` are Wasserstein distances between those transformed values.
#'
#' A positive `J2` means higher average squared nearest-observation correlation
#' during validation. Under Gaussian prediction using one neighbour with known
#' mean and covariance, target MSE minus expected CV MSE equals total variance
#' times `J2`. This is not a general bias correction for fitted prediction models.
#' The integral interpretation of `U1` and `U2` as unsigned weighted distance-CDF
#' contrasts assumes nonnegative correlation decreasing with distance.
#'
#' @param validation,target Nonempty numeric vectors of finite, nonnegative
#'   nearest-training distances.
#' @param correlation A function mapping distances to correlations in `[-1, 1]`,
#'   or a [spatial_covariance()] object. If `NULL`, correlation diagnostics are `NA`.
#' @param validation_weights,target_weights Nonnegative weights, normalised
#'   separately. Defaults are equal weights within each distribution.
#' @return A named numeric vector: `A`, `W`, `J1`, `J2`, `U1`, `U2`, `Aplus`,
#'   `Aminus`, and `crosses` (one if both positive and negative CDF areas exceed
#'   a relative floating-point tolerance, otherwise zero).
#' @export
#' @examples
#' model <- spatial_covariance(range = 1, signal = 0.8)
#' distance_contrasts(c(0.2, 0.5, 1), c(0.5, 1, 2), model)
distance_contrasts <- function(validation, target, correlation = NULL,
                               validation_weights = NULL, target_weights = NULL) {
  for (nm in c("validation", "target")) {
    z <- get(nm)
    if (!is.numeric(z) || !length(z) || any(!is.finite(z)) || any(z < 0))
      stop("`", nm, "` must contain finite nonnegative distances.", call. = FALSE)
  }
  wc <- .weights(validation_weights, length(validation), "validation_weights")
  wq <- .weights(target_weights, length(target), "target_weights")
  cfun <- .correlation(correlation)
  ar <- .cdf_areas(validation, target, wc, wq)
  j <- c(J1 = NA_real_, J2 = NA_real_, U1 = NA_real_, U2 = NA_real_)
  if (!is.null(cfun)) {
    rc <- cfun(validation); rq <- cfun(target)
    if (!is.numeric(rc) || !is.numeric(rq) || length(rc) != length(validation) ||
        length(rq) != length(target) || any(!is.finite(c(rc, rq))) || any(abs(c(rc, rq)) > 1))
      stop("The correlation function must return one finite value in [-1, 1] per distance.", call. = FALSE)
    j <- c(J1 = sum(wc * rc) - sum(wq * rq), J2 = sum(wc * rc^2) - sum(wq * rq^2),
           U1 = unname(.cdf_areas(rc, rq, wc, wq)["W"]),
           U2 = unname(.cdf_areas(rc^2, rq^2, wc, wq)["W"]))
  }
  c(A = sum(wq * target) - sum(wc * validation), W = unname(ar["W"]), j,
    Aplus = unname(ar["positive"]), Aminus = unname(ar["negative"]),
    crosses = as.numeric(all(ar[1:2] > 64 * .Machine$double.eps * max(ar["W"], .Machine$double.eps))))
}
