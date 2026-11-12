.scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper || (integer && x != floor(x)))
    stop("`", name, "` must be a finite ", if (integer) "integer" else "number",
         " between ", lower, " and ", upper, ".", call. = FALSE)
  invisible(x)
}

.weights <- function(w, n, name = "weights") {
  if (is.null(w)) return(rep(1 / n, n))
  if (!is.numeric(w) || length(w) != n || any(!is.finite(w)) ||
      any(w < 0) || max(w) <= 0)
    stop("`", name, "` must have one finite nonnegative value per observation and a positive sum.",
         call. = FALSE)
  w <- w / max(w)
  w / sum(w)
}

.coordinates <- function(x, name = "locations") {
  crs <- NULL
  if (inherits(x, "sf") || inherits(x, "sfc")) {
    if (!requireNamespace("sf", quietly = TRUE))
      stop("Install 'sf' to use spatial objects.", call. = FALSE)
    if (is.na(sf::st_crs(x)) || isTRUE(sf::st_is_longlat(x)))
      stop("`", name, "` must have a projected CRS. Transform longitude/latitude with sf::st_transform().",
           call. = FALSE)
    if (any(sf::st_geometry_type(x) != "POINT") || any(sf::st_is_empty(x)))
      stop("`", name, "` must contain nonempty POINT geometries.", call. = FALSE)
    crs <- sf::st_crs(x)
    x <- sf::st_coordinates(x)[, 1:2, drop = FALSE]
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) != 2L ||
      nrow(x) < 1L || any(!is.finite(x)))
    stop("`", name, "` must be a finite numeric matrix with two coordinate columns, or projected sf points.",
         call. = FALSE)
  list(xy = unname(x), crs = crs)
}

.target_coordinates <- function(locations, target) {
  x <- .coordinates(locations)
  if (!inherits(target, "spcv_target")) target <- prediction_target(target)
  if (xor(is.null(x$crs), is.null(target$crs)))
    stop("Use either numeric coordinates for both locations and target, or projected sf objects for both.",
         call. = FALSE)
  if (!is.null(x$crs) && !isTRUE(x$crs == target$crs))
    stop("Training and target CRS differ; transform them to the same projected CRS first.", call. = FALSE)
  list(x = x$xy, q = target$locations, target = target, crs = x$crs)
}

.pair_distance <- function(x, y) {
  # Differences avoid catastrophic cancellation for large projected coordinates.
  sqrt(outer(x[, 1L], y[, 1L], "-")^2 + outer(x[, 2L], y[, 2L], "-")^2)
}

.nearest <- function(x, y) {
  out <- numeric(nrow(x))
  for (start in seq.int(1L, nrow(x), by = 256L)) {
    ii <- seq.int(start, min(nrow(x), start + 255L))
    out[ii] <- apply(.pair_distance(x[ii, , drop = FALSE], y), 1L, min)
  }
  out
}

.local_seed <- function(seed, expr) {
  if (is.null(seed)) return(force(expr))
  .scalar(seed, "seed", 0, .Machine$integer.max, TRUE)
  exists_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (exists_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (exists_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  force(expr)
}

.correlation <- function(covariance) {
  if (inherits(covariance, "spcv_covariance")) return(function(h) spatial_correlation(covariance, h))
  if (is.function(covariance)) return(covariance)
  if (is.null(covariance)) return(NULL)
  stop("`covariance` must be a spatial_covariance() object or a correlation function.", call. = FALSE)
}
