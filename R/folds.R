#' Validate an explicit cross-validation partition
#'
#' Every observation must be tested exactly once. Training and test indices must
#' be nonempty, unique, disjoint, and within `1:n`. Training sets may exclude
#' additional observations, as in buffered CV. Repeated CV should be audited one
#' complete partition at a time.
#'
#' @param folds A list whose elements each contain integer vectors `train` and
#'   `test`; or a vector of fold labels of length `n`.
#' @param n Number of training observations.
#' @return A list of class `spcv_folds`, with one `train`/`test` list per fold and
#'   an `n` attribute recording the number of observations.
#' @export
#' @examples
#' cv_folds(c(1, 2, 1, 2), n = 4)
cv_folds <- function(folds, n) {
  .scalar(n, "n", 2, .Machine$integer.max, TRUE)
  if (!is.list(folds)) {
    if (length(folds) != n || anyNA(folds) ||
        !(is.numeric(folds) || is.character(folds) || is.factor(folds)))
      stop("Fold labels must have length `n` and no missing values.", call. = FALSE)
    folds <- lapply(unique(folds), function(id)
      list(train = which(folds != id), test = which(folds == id)))
  }
  if (!length(folds)) stop("At least one fold is required.", call. = FALSE)
  index <- function(z) is.numeric(z) && length(z) > 0 && all(is.finite(z)) &&
    all(z == floor(z)) && all(z >= 1 & z <= n) && !anyDuplicated(z)
  out <- lapply(seq_along(folds), function(i) {
    f <- folds[[i]]
    if (!is.list(f) || !index(f$train) || !index(f$test))
      stop("Fold ", i, " must contain nonempty, unique integer `train` and `test` indices within 1:n.", call. = FALSE)
    if (length(intersect(f$train, f$test))) stop("Fold ", i, " has overlapping training and test indices.", call. = FALSE)
    list(train = as.integer(f$train), test = as.integer(f$test))
  })
  tested <- unlist(lapply(out, `[[`, "test"), use.names = FALSE)
  if (length(tested) != n || any(tabulate(tested, nbins = n) != 1L))
    stop("Every observation must occur in exactly one test fold.", call. = FALSE)
  structure(out, class = c("spcv_folds", "list"), n = as.integer(n))
}

#' Construct random, spatial-block, buffered, or kNNDM folds
#'
#' Spatial blocks use coordinate k-means. Buffered folds remove training
#' observations within `buffer` of any test observation. The kNNDM method calls
#' [CAST::knndm()] and requires the optional packages 'CAST' and 'sf'. Its target
#' must have equal weights. All methods return explicit indices suitable for
#' [audit_cv()]. These constructors do not use response values.
#'
#' @param locations Numeric two-column coordinates or projected spatial points.
#' @param k Number of folds, between two and the number of observations.
#' @param method One of `"random"`, `"block"`, `"buffered"`, or `"knndm"`.
#' @param buffer Nonnegative exclusion distance, required for buffered CV and
#'   expressed in the same units as the coordinates.
#' @param target Prediction target, required for kNNDM.
#' @param seed An integer seed. The default is reproducible and restores the
#'   caller's random-number state, including after calls to external constructors.
#'   Use `NULL` to use and advance the caller's random-number stream.
#' @return Validated [cv_folds()] with explicit training and test indices.
#' @export
#' @examples
#' xy <- cbind(seq_len(12), rep(c(0, 1), 6))
#' spatial_folds(xy, k = 3, method = "random", seed = 10)
spatial_folds <- function(locations, k = 5L,
                          method = c("random", "block", "buffered", "knndm"),
                          buffer = NULL, target = NULL, seed = 1L) {
  x <- .coordinates(locations)$xy
  n <- nrow(x)
  .scalar(k, "k", 2, n, TRUE)
  method <- match.arg(method)
  if (method == "buffered") .scalar(buffer, "buffer", 0)
  .local_seed(seed, {
    if (method == "knndm") {
      if (!requireNamespace("CAST", quietly = TRUE) || !requireNamespace("sf", quietly = TRUE))
        stop("Install 'CAST' and 'sf' to construct kNNDM folds, or supply external folds to cv_folds().", call. = FALSE)
      if (is.null(target)) stop("kNNDM requires a prediction target.", call. = FALSE)
      z <- .target_coordinates(locations, target)
      if (max(z$target$weights) - min(z$target$weights) > 64 * .Machine$double.eps)
        stop("The CAST kNNDM constructor requires equal target weights; weighted targets can still be audited with audit_cv().", call. = FALSE)
      # CAST needs a projected CRS to request Euclidean distances. Numeric input
      # is assigned a planar CRS solely for calculation, without changing units.
      crs <- if (is.null(z$crs)) 3857 else z$crs
      points <- function(m) sf::st_as_sf(as.data.frame(m), coords = 1:2, crs = crs)
      fit <- CAST::knndm(points(x), predpoints = points(z$q), k = k)
      result <- Map(function(tr, te) list(train = tr, test = te), fit$indx_train, fit$indx_test)
    } else {
      ids <- if (method == "random") sample(rep_len(seq_len(k), n)) else {
        if (nrow(unique(x)) < k) stop("Spatial blocks need at least k distinct locations.", call. = FALSE)
        if (k == n) seq_len(n) else stats::kmeans(x, centers = k, nstart = 10L, iter.max = 100L)$cluster
      }
      result <- cv_folds(ids, n)
      if (method == "buffered") result <- lapply(result, function(f) {
        d <- .nearest(x[f$train, , drop = FALSE], x[f$test, , drop = FALSE])
        f$train <- f$train[d > buffer]
        f
      })
    }
    cv_folds(result, n)
  })
}
