#' Audit a CV partition against a prediction target
#'
#' Each validation observation is compared with its fold's available training
#' observations. Each target location is compared with the full training set.
#' Validation observations receive weight `1/n`; target weights come from
#' [prediction_target()]. This function needs no response values.
#'
#' @param locations Numeric two-column coordinates or projected spatial points.
#' @param folds Fold labels or an explicit list accepted by [cv_folds()].
#' @param target A [prediction_target()] object, or target coordinates.
#' @param covariance A [spatial_covariance()] or [estimate_covariance()] object,
#'   a correlation function, or `NULL` for distance diagnostics only.
#' @return An object of class `spcv_audit` containing `diagnostics`, validation
#'   and target nearest distances, normalised target weights, validated folds,
#'   and the covariance specification. `diagnostics` is described in
#'   [distance_contrasts()].
#' @export
#' @examples
#' xy <- cbind(seq_len(12) / 12, rep(c(0, 1), 6))
#' q <- prediction_target(rbind(c(0.2, 0.5), c(0.8, 0.5)))
#' folds <- spatial_folds(xy, k = 3)
#' audit_cv(xy, folds, q, spatial_covariance(0.3, 0.8))
audit_cv <- function(locations, folds, target, covariance = NULL) {
  z <- .target_coordinates(locations, target)
  folds <- cv_folds(folds, nrow(z$x))
  dc <- numeric(nrow(z$x))
  for (f in folds)
    dc[f$test] <- .nearest(z$x[f$test, , drop = FALSE], z$x[f$train, , drop = FALSE])
  dq <- .nearest(z$q, z$x)
  diagnostics <- distance_contrasts(dc, dq, covariance, target_weights = z$target$weights)
  structure(list(diagnostics = diagnostics, validation_distance = dc,
                 target_distance = dq, target_weights = z$target$weights,
                 folds = folds, covariance = covariance, target = z$target),
            class = "spcv_audit")
}

#' Compare candidate CV schemes using a common diagnostic
#'
#' `A`, `J1` and `J2` are ranked by absolute value; `W`, `U1` and `U2` by value.
#' All candidates use the same supplied covariance and target. Schemes within
#' `tolerance` of the minimum are tied and receive equal selection probabilities.
#' No observed target errors are used to select or break ties. When evaluating
#' such a rule, average the individual tied schemes' losses with these
#' probabilities; do not average their signed errors before taking absolute values.
#'
#' @param locations Numeric two-column coordinates or projected spatial points.
#' @param schemes A nonempty named list of candidate fold partitions.
#' @param target A [prediction_target()] or target coordinates.
#' @param covariance A covariance object or correlation function; required for
#'   correlation-based selection.
#' @param criterion One of `"J2"`, `"W"`, `"J1"`, `"A"`, `"U1"`, or `"U2"`.
#' @param tolerance Absolute numerical tolerance for ties, in criterion units.
#' @return A `spcv_comparison` list with `table` (one row per scheme, diagnostics,
#'   `criterion_value`, `selected`, and `selection_probability`), individual
#'   `audits`, the selection criterion and tolerance.
#' @export
#' @examples
#' xy <- cbind(seq_len(20) / 20, rep(c(0, 1), 10))
#' q <- prediction_target(rbind(c(0.2, 0.5), c(0.8, 0.5)))
#' schemes <- list(random = spatial_folds(xy, 4),
#'                 block = spatial_folds(xy, 4, "block"))
#' compare_cv_schemes(xy, schemes, q, spatial_covariance(0.3, 0.8))
compare_cv_schemes <- function(locations, schemes, target, covariance = NULL,
                               criterion = c("J2", "W", "J1", "A", "U1", "U2"),
                               tolerance = 1e-12) {
  criterion <- match.arg(criterion)
  .scalar(tolerance, "tolerance", 0)
  if (!is.list(schemes) || !length(schemes) || is.null(names(schemes)) ||
      anyNA(names(schemes)) || any(names(schemes) == "") || anyDuplicated(names(schemes)))
    stop("`schemes` must be a nonempty list with unique, nonempty names.", call. = FALSE)
  audits <- lapply(schemes, function(f) audit_cv(locations, f, target, covariance))
  d <- as.data.frame(do.call(rbind, lapply(audits, `[[`, "diagnostics")))
  value <- if (criterion %in% c("A", "J1", "J2")) abs(d[[criterion]]) else d[[criterion]]
  if (anyNA(value)) stop("Supply a covariance object or correlation function for this criterion.", call. = FALSE)
  selected <- value <= min(value) + tolerance
  tab <- data.frame(scheme = names(schemes), d, criterion_value = value,
                    selected = selected, selection_probability = selected / sum(selected),
                    row.names = NULL, check.names = FALSE)
  structure(list(table = tab, audits = audits, criterion = criterion,
                 tolerance = tolerance), class = "spcv_comparison")
}

#' @export
print.spcv_audit <- function(x, ...) {
  cat("Spatial CV audit\n")
  print(x$diagnostics[c("A", "W", "J1", "J2", "U1", "U2")])
  cat("Diagnostics compare prediction conditions; they do not certify CV accuracy.\n")
  invisible(x)
}

#' @export
summary.spcv_audit <- function(object, ...) {
  c(object$diagnostics, n_validation = length(object$validation_distance),
    n_target = length(object$target_distance), n_folds = length(object$folds))
}

#' @export
print.spcv_comparison <- function(x, ...) {
  cat("CV comparison by", if (x$criterion %in% c("A", "J1", "J2"))
    paste0("|", x$criterion, "|") else x$criterion, "\n")
  print(x$table[, c("scheme", "W", "J2", "criterion_value", "selection_probability")], row.names = FALSE)
  invisible(x)
}

#' Plot validation and target distance distributions
#'
#' The solid curve is validation and the dashed curve is the weighted target.
#' @param x A [audit_cv()] result.
#' @param ... Additional arguments passed to [graphics::plot()].
#' @return Invisibly, the audit object. Called for its plotting side effect.
#' @export
#' @examples
#' xy <- cbind(1:8, rep(c(0, 1), 4))
#' a <- audit_cv(xy, rep(1:2, 4), rbind(c(2, 2), c(6, 2)))
#' plot(a)
plot.spcv_audit <- function(x, ...) {
  cuts <- sort(unique(c(0, x$validation_distance, x$target_distance)))
  fc <- stats::ecdf(x$validation_distance)(cuts)
  fq <- vapply(cuts, function(d) sum(x$target_weights[x$target_distance <= d]), numeric(1))
  args <- list(x = cuts, y = fc, type = "s", ylim = c(0, 1),
               xlab = "Nearest training distance", ylab = "Cumulative proportion",
               col = "#146C75", lwd = 2, lty = 1)
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(graphics::plot, args)
  graphics::lines(cuts, fq, type = "s", col = "#C46530", lwd = 2, lty = 2)
  graphics::legend("bottomright", c("Validation (solid)", "Target (dashed)"),
                   col = c("#146C75", "#C46530"), lty = c(1, 2), lwd = 2, bty = "n")
  invisible(x)
}
