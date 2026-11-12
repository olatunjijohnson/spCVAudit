test_that("weighted empirical contrasts agree with hand calculations", {
  r <- distance_contrasts(c(0, 2), c(1, 3), function(d) exp(-d),
                          validation_weights = c(1, 3), target_weights = c(3, 1))
  expect_equal(unname(r["A"]), 0)
  expect_equal(unname(r["W"]), 1)
  expect_equal(unname(r["Aplus"]), 0.5)
  expect_equal(unname(r["Aminus"]), 0.5)
  expect_equal(unname(r["J2"]), (1 + 3 * exp(-4) - 3 * exp(-2) - exp(-6)) / 4)
  expect_equal(unname(r["crosses"]), 1)
  expect_lte(abs(r[["J1"]]), r[["U1"]] + 1e-15)
  expect_lte(abs(r[["J2"]]), r[["U2"]] + 1e-15)
})

test_that("distance integration is exact at atoms and respects weights", {
  expect_equal(distance_contrasts(c(1, 1), 1)[["W"]], 0)
  expect_equal(distance_contrasts(c(0, 100), 2, validation_weights = c(1, 0))[["W"]], 2)
  a <- distance_contrasts(c(1, 2), c(2, 5), spatial_covariance(2, 0.4))
  b <- distance_contrasts(c(2, 5), c(1, 2), spatial_covariance(2, 0.4))
  expect_equal(a[c("A", "J1", "J2")], -b[c("A", "J1", "J2")])
  expect_equal(a[c("W", "U1", "U2")], b[c("W", "U1", "U2")])
  c <- distance_contrasts(c(1, 2), c(2, 5), spatial_covariance(2, 0.8))
  expect_equal(c[["J2"]], 4 * a[["J2"]])
  expect_equal(c[["J1"]], 2 * a[["J1"]])
  expect_error(distance_contrasts(-1, 2), "nonnegative")
  expect_error(distance_contrasts(1, 2, function(d) 2), "correlation")
  expect_error(distance_contrasts(1, 2, target_weights = 0), "positive sum")
  expect_error(distance_contrasts(numeric(), 2), "distances")
})

test_that("J1 and J2 can have opposite signs", {
  d <- distance_contrasts(-log(0.5), -log(c(0.8, 0.01)),
                          spatial_covariance(1), target_weights = c(0.6, 0.4))
  expect_equal(d[["J1"]], 0.016)
  expect_equal(d[["J2"]], -0.13404)
})

test_that("invalid targets and covariance parameters are rejected", {
  expect_error(prediction_target(matrix(1, 2, 3)), "two coordinate")
  expect_error(prediction_target(matrix(c(NA, 1), 1)), "finite")
  expect_error(prediction_target(matrix(1:4, 2), c(-1, 2)), "nonnegative")
  expect_error(spatial_covariance(0), "range")
  expect_error(spatial_covariance(1, 1.1), "signal")
  expect_error(spatial_covariance(1, variance = 0), "variance")
  expect_equal(spatial_correlation(spatial_covariance(1, 0.8), 0), 0.8)
  expect_equal(spatial_correlation(spatial_covariance(1e-300, family = "matern32"), 1e300), 0)
})
