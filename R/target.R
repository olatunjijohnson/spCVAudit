#' Declare the locations and weights of a prediction target
#'
#' The target describes where prediction error is to be assessed. Target response
#' values are not required and never enter the diagnostics. Coordinates use a
#' common projected system; numeric matrices are assumed to share their units.
#'
#' @param locations A numeric matrix or data frame with exactly two coordinate
#'   columns, or an `sf`/`sfc` object of points in a projected CRS.
#' @param weights Nonnegative target weights, normalised to sum to one. The
#'   default gives each target location equal weight.
#' @param label A single descriptive character string.
#' @return An object of class `spcv_target` containing locations, normalised
#'   weights, label and, for spatial objects, CRS.
#' @export
#' @examples
#' target <- prediction_target(rbind(c(0, 0), c(1, 1)), weights = c(1, 3))
#' target
prediction_target <- function(locations, weights = NULL, label = "Prediction target") {
  x <- .coordinates(locations)
  if (!is.character(label) || length(label) != 1L || is.na(label))
    stop("`label` must be one character string.", call. = FALSE)
  structure(list(locations = x$xy, weights = .weights(weights, nrow(x$xy)),
                 label = label, crs = x$crs), class = "spcv_target")
}

#' @export
print.spcv_target <- function(x, ...) {
  cat(x$label, "\n", nrow(x$locations), " locations; ",
      if (length(unique(x$weights)) == 1L) "equal" else "unequal",
      " weights\n", sep = "")
  invisible(x)
}
