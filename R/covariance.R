#' Specify a stationary isotropic spatial covariance
#'
#' `range` is the scale parameter, not a common practical correlation range across
#' families. For distinct observations, correlation is `signal * exp(-d/range)`
#' for the exponential family and `signal * (1 + sqrt(3)*d/range) *
#' exp(-sqrt(3)*d/range)` for the Matern-3/2 family. Independent observation noise
#' has variance `variance * (1 - signal)`, including at coincident locations.
#'
#' @param range Positive correlation scale in coordinate-distance units.
#' @param signal Spatial signal fraction in `[0, 1]`.
#' @param variance Positive total response variance.
#' @param family Correlation family, `"exponential"` or `"matern32"`.
#' @return A `spcv_covariance` object with the supplied parameters, spatial
#'   variance `psill`, independent noise variance `nugget`, and `estimated = FALSE`.
#' @export
#' @examples
#' model <- spatial_covariance(range = 0.2, signal = 0.8)
#' spatial_correlation(model, c(0, 0.2, 1))
spatial_covariance <- function(range, signal = 1, variance = 1,
                               family = c("exponential", "matern32")) {
  .scalar(range, "range", .Machine$double.xmin)
  .scalar(signal, "signal", 0, 1)
  .scalar(variance, "variance", .Machine$double.xmin)
  structure(list(range = range, signal = signal, variance = variance,
                 family = match.arg(family), psill = variance * signal,
                 nugget = variance * (1 - signal), estimated = FALSE),
            class = "spcv_covariance")
}

#' Evaluate correlation between distinct spatial observations
#' @param model A [spatial_covariance()] or [estimate_covariance()] object.
#' @param distance Finite, nonnegative distances.
#' @return A numeric vector with one correlation per distance. At zero distance
#'   this is the signal fraction, because the observations have independent noise.
#' @export
#' @examples
#' spatial_correlation(spatial_covariance(1, 0.8), c(0, 1, 2))
spatial_correlation <- function(model, distance) {
  if (!inherits(model, "spcv_covariance")) stop("`model` must be a spatial covariance object.", call. = FALSE)
  if (!is.numeric(distance) || any(!is.finite(distance)) || any(distance < 0))
    stop("`distance` must be finite and nonnegative.", call. = FALSE)
  z <- distance / model$range
  r <- if (model$family == "exponential") exp(-z) else {
    z <- sqrt(3) * z
    out <- (1 + z) * exp(-z)
    out[is.infinite(z)] <- 0
    out
  }
  model$signal * r
}

#' @export
print.spcv_covariance <- function(x, ...) {
  cat(if (isTRUE(x$estimated)) "Estimated" else "Specified", x$family, "covariance\n")
  cat("  scale:", format(x$range), " signal:", format(x$signal),
      " total variance:", format(x$variance), "\n")
  if (isTRUE(x$estimated)) cat("  fit boundary:", paste(x$fit$boundary, collapse = ", "), "\n")
  invisible(x)
}
