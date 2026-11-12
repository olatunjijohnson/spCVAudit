test_that("profile fit is reproducible, scale equivariant and feasible", {
  xy <- as.matrix(expand.grid(x = seq(0, 1, length.out = 7), y = seq(0, 1, length.out = 7)))
  y <- sin(3 * xy[, 1]) + cos(5 * xy[, 2])
  set.seed(7)
  before <- .Random.seed
  fit <- estimate_covariance(xy, y)
  expect_identical(.Random.seed, before)
  expect_true(fit$estimated)
  expect_gte(fit$signal, 0)
  expect_lte(fit$signal, 1)
  expect_gt(fit$variance, 0)
  scaled <- estimate_covariance(xy * 1000, y * 3)
  expect_equal(scaled$range / 1000, fit$range, tolerance = 1e-5)
  expect_equal(scaled$signal, fit$signal, tolerance = 1e-5)
  expect_equal(scaled$variance / 9, fit$variance, tolerance = 1e-5)
  # Profile objective cannot exceed any feasible tested parameter combination.
  vg <- fit$variogram
  for (range in c(0.01, 0.1, 1)) {
    trial <- 0.1 + 0.5 * (1 - exp(-vg$distance / range))
    loss <- weighted.mean((trial - vg$semivariance)^2, vg$pairs)
    expect_lte(fit$fit$objective, loss + 1e-10)
  }
  expect_error(estimate_covariance(xy, rep(1, nrow(xy))), "constant")
  expect_error(estimate_covariance(xy, c(NA, y[-1])), "finite")
  expect_error(estimate_covariance(xy[1:4, ], y[1:4]), "five")
})

test_that("bounded pair sampling produces distinct valid unordered pairs", {
  set.seed(42)
  p <- spCVAudit:::.sample_pairs(1000, 200)
  expect_true(all(p[1, ] < p[2, ]))
  expect_true(all(p >= 1 & p <= 1000))
  expect_equal(anyDuplicated(t(p)), 0L)
})
