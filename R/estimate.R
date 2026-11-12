.sample_pairs <- function(n, max_pairs) {
  total <- as.double(n) * (n - 1) / 2
  if (total <= max_pairs) return(utils::combn(n, 2L))
  ranks <- sample.int(total, max_pairs)
  i <- floor(n - 0.5 - sqrt(pmax((n - 0.5)^2 - 2 * (ranks - 1), 0))) + 1
  i <- pmin(pmax(i, 1), n - 1)
  start <- (i - 1) * (2 * n - i) / 2
  i[ranks <= start] <- i[ranks <= start] - 1
  start <- (i - 1) * (2 * n - i) / 2
  i[ranks > start + n - i] <- i[ranks > start + n - i] + 1
  start <- (i - 1) * (2 * n - i) / 2
  rbind(i, ranks - start + i)
}

#' Estimate exponential covariance by profiling a binned variogram
#'
#' Fits `nugget + psill * (1 - exp(-distance / range))` to an empirical
#' semivariogram, with bin weights proportional to pair counts. For each range,
#' constrained nugget and partial sill are obtained by evaluating the quadratic
#' optimum and all parameter boundaries. An 81-point log-range grid and local
#' refinements reduce dependence on optimiser starting values. This does not
#' account for covariance estimation uncertainty or establish model adequacy.
#'
#' The fit is isotropic. Supply residuals as `response` to describe residual
#' rather than marginal dependence. To avoid leakage in CV predictions, estimate
#' covariance separately using each fold's training responses. A common estimate
#' from all training responses may be used to compare diagnostics across schemes;
#' that selection then depends on responses as well as locations.
#'
#' @param locations Numeric two-column coordinates or projected spatial points.
#' @param response Finite numeric observations or residuals, one per location.
#' @param n_bins Number of equal-width distance bins.
#' @param cutoff Maximum pair distance, default half the maximum sampled distance.
#' @param max_pairs Maximum number of unordered pairs; larger datasets use a
#'   sample without replacement, before allocating pair coordinates.
#' @param seed Seed for pair sampling; restores the caller's random-number state.
#' @return A `spcv_covariance` object with `estimated = TRUE`, its empirical
#'   `variogram`, and `fit` diagnostics including objective, boundary flags and
#'   scaled parameter limits. Boundary fits remain visible for assessment.
#' @export
#' @examples
#' xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = 6),
#'                            y = seq(0, 1, length.out = 6)))
#' y <- sin(3 * xy[, 1]) + cos(5 * xy[, 2])
#' fit <- estimate_covariance(xy, y)
#' fit
estimate_covariance <- function(locations, response, n_bins = 15L, cutoff = NULL,
                                 max_pairs = 100000L, seed = 1L) {
  xy <- .coordinates(locations)$xy
  n <- nrow(xy)
  if (n < 5L) stop("At least five observations are required to estimate covariance.", call. = FALSE)
  if (!is.numeric(response) || length(response) != n || any(!is.finite(response)))
    stop("`response` must contain one finite numeric value per location.", call. = FALSE)
  .scalar(n_bins, "n_bins", 3, 10000, TRUE)
  .scalar(max_pairs, "max_pairs", 3, .Machine$integer.max, TRUE)
  pairs <- .local_seed(seed, .sample_pairs(n, max_pairs))
  d <- sqrt(rowSums((xy[pairs[1, ], , drop = FALSE] - xy[pairs[2, ], , drop = FALSE])^2))
  gamma <- 0.5 * (response[pairs[1, ]] - response[pairs[2, ]])^2
  if (any(!is.finite(gamma))) stop("Response differences are too large for stable variogram calculation.", call. = FALSE)
  if (is.null(cutoff)) cutoff <- max(d) / 2
  .scalar(cutoff, "cutoff", .Machine$double.xmin)
  # Work on relative distances so changing coordinate units cannot move a pair
  # across a bin edge through last-bit rounding (especially on regular grids).
  relative <- round(d / cutoff, 12L)
  keep <- d > 0 & relative <= 1
  bin <- factor(pmax(1L, pmin(n_bins, ceiling(round(relative[keep] * n_bins, 10L)))),
                levels = seq_len(n_bins))
  h <- tapply(d[keep], bin, mean)
  g <- tapply(gamma[keep], bin, mean)
  np <- tapply(gamma[keep], bin, length)
  ok <- !is.na(h)
  vg <- data.frame(distance = as.numeric(h[ok]), semivariance = as.numeric(g[ok]),
                    pairs = as.integer(np[ok]))
  if (nrow(vg) < 3L) stop("At least three occupied distance bins are needed; adjust n_bins or cutoff.", call. = FALSE)
  hs <- max(vg$distance); gs <- max(vg$semivariance)
  if (gs <= 0) stop("Cannot estimate a positive covariance from constant responses.", call. = FALSE)
  h <- vg$distance / hs; y <- vg$semivariance / gs; w <- vg$pairs / sum(vg$pairs)
  ey <- sum(w * y)
  clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)
  solve_range <- function(log_range, coefficients = FALSE) {
    x <- -expm1(-h / exp(log_range))
    ex <- sum(w * x); exx <- sum(w * x^2); exy <- sum(w * x * y)
    candidates <- rbind(c(0, clip(exy / exx, 0, 2)),
                        c(1, clip((exy - ex) / exx, 0, 2)),
                        c(clip(ey, 0, 1), 0), c(clip(ey - 2 * ex, 0, 1), 2))
    if (exx - ex^2 > 1e-15) {
      b <- (exy - ex * ey) / (exx - ex^2); a <- ey - b * ex
      if (a >= 0 && a <= 1 && b >= 0 && b <= 2) candidates <- rbind(candidates, c(a, b))
    }
    loss <- apply(candidates, 1L, function(p) sum(w * (p[1] + p[2] * x - y)^2))
    best <- which.min(loss)
    if (!coefficients) return(loss[best])
    c(nugget = candidates[best, 1], psill = candidates[best, 2],
      range = exp(log_range), objective = loss[best])
  }
  grid <- seq(log(1e-6), log(5), length.out = 81L)
  loss <- vapply(grid, solve_range, numeric(1))
  local <- which(loss[2:80] <= loss[1:79] & loss[2:80] <= loss[3:81]) + 1L
  local <- local[vapply(local, function(i) loss[i] < max(loss[i - 1L], loss[i + 1L]) - 1e-14, logical(1))]
  trials <- c(grid[which.min(loss)], grid[1], tail(grid, 1))
  for (i in local)
    trials <- c(trials, stats::optimize(solve_range, grid[c(i - 1L, i + 1L)], tol = 1e-10)$minimum)
  solutions <- lapply(trials, solve_range, coefficients = TRUE)
  best <- solutions[[which.min(vapply(solutions, `[[`, numeric(1), "objective"))]]
  if (sum(best[c("nugget", "psill")]) <= 0) stop("The fitted total variance is zero.", call. = FALSE)
  out <- spatial_covariance(unname(best["range"] * hs),
                            unname(best["psill"] / sum(best[c("nugget", "psill")])),
                            unname(sum(best[c("nugget", "psill")]) * gs))
  boundary <- c(if (best["nugget"] <= 1e-8) "zero nugget",
                if (best["psill"] <= 1e-8) "zero spatial signal",
                if (best["nugget"] >= 1 - 1e-8) "nugget upper limit",
                if (best["psill"] >= 2 - 1e-8) "spatial variance upper limit",
                if (best["range"] <= 1e-6 * (1 + 1e-6)) "range lower limit",
                if (best["range"] >= 5 * (1 - 1e-6)) "range upper limit")
  out$estimated <- TRUE
  out$variogram <- vg
  out$fit <- list(objective = unname(best["objective"] * gs^2),
                   boundary = if (length(boundary)) boundary else "none",
                   method = "Profile weighted least squares", seed = seed,
                   pairs_used = sum(vg$pairs), cutoff = cutoff,
                   limits = list(range = c(1e-6, 5) * hs, nugget = c(0, 1) * gs,
                                 psill = c(0, 2) * gs))
  out
}
