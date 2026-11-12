test_that("paper counterexample has zero nearest contrast and nonzero full bias", {
  xy <- cbind(c(-1, 0, 1), 0)
  q <- matrix(c(2, 0), 1)
  r <- gaussian_cv(xy, 1:3, q, spatial_covariance(1, 0.8))
  expect_equal(r$J2, 0)
  expect_lte(abs(r$optimism - 0.02305608), 5e-9)
  expect_equal(r$optimism, r$one_neighbour + r$additional_neighbours)
  expect_lte(r$lower, r$optimism)
  expect_gte(r$upper, r$optimism)
})

test_that("full MSE agrees with independent direct matrix calculations", {
  set.seed(51)
  for (family in c("exponential", "matern32")) {
    x <- matrix(runif(16), 8, 2)
    q <- matrix(runif(6), 3, 2)
    w <- c(0.2, 0.3, 0.5)
    model <- spatial_covariance(0.4, 0.7, 2.5, family)
    rho <- function(d) if (family == "exponential") exp(-d / 0.4) else
      (1 + sqrt(3) * d / 0.4) * exp(-sqrt(3) * d / 0.4)
    direct <- function(train, query) {
      d <- as.matrix(stats::dist(rbind(train, query)))
      n <- nrow(train)
      s <- 1.75 * rho(d[seq_len(n), seq_len(n), drop = FALSE]) + diag(0.75, n)
      b <- 1.75 * rho(d[n + seq_len(nrow(query)), seq_len(n), drop = FALSE])
      2.5 - diag(b %*% solve(s, t(b)))
    }
    folds <- cv_folds(rep(1:4, 2), 8)
    mse <- numeric(8)
    for (f in folds) mse[f$test] <- direct(x[f$train, , drop = FALSE], x[f$test, , drop = FALSE])
    a <- gaussian_cv(x, folds, prediction_target(q, w), model)
    expect_equal(a$cv_mse, mean(mse), tolerance = 1e-12)
    expect_equal(a$target_mse, sum(w * direct(x, q)), tolerance = 1e-12)
    expect_equal(a$optimism, a$one_neighbour + a$additional_neighbours, tolerance = 1e-12)
    expect_equal(a$J2, audit_cv(x, folds, prediction_target(q, w), model)$diagnostics[["J2"]])
    expect_true(all(a$validation$additional <= a$validation$bound + 1e-12))
    expect_true(all(a$target$additional <= a$target$bound + 1e-12))
  }
})

test_that("bound holds across signal fractions, locations and families", {
  set.seed(984)
  for (i in seq_len(40)) {
    x <- matrix(runif(12), 6, 2)
    q <- matrix(runif(8), 4, 2)
    m <- spatial_covariance(runif(1, 0.02, 2), runif(1, 0, 0.999),
                            runif(1, 0.1, 5), sample(c("exponential", "matern32"), 1))
    r <- gaussian_cv(x, rep(1:3, 2), prediction_target(q, runif(4)), m)
    expect_lte(r$lower, r$optimism + 1e-10)
    expect_gte(r$upper, r$optimism - 1e-10)
  }
})

test_that("independent noise and coincident locations have correct semantics", {
  x <- cbind(c(0, 0, 1), 0)
  q <- matrix(c(0, 0), 1)
  r <- gaussian_cv(x, 1:3, q, spatial_covariance(1, 0, variance = 2))
  expect_equal(r$target_mse, 2)
  expect_equal(r$cv_mse, 2)
  expect_equal(r$optimism, 0)
  noisy <- gaussian_cv(x, 1:3, q, spatial_covariance(1, 0.8))
  expect_gte(noisy$target_mse, 0.2 - 1e-12)
  expect_error(gaussian_cv(x, 1:3, q, spatial_covariance(1, 1)), "positive definite")
  zero_nugget <- gaussian_cv(cbind(1:3, 0), 1:3, q, spatial_covariance(1, 1))
  expect_true(is.na(zero_nugget$lower))
})

test_that("conditional target MSE includes mean discrepancy", {
  x <- cbind(1:3, 0)
  q <- matrix(c(1.5, 0), 1)
  m <- spatial_covariance(1, 0.8)
  r <- gaussian_target_mse(x, c(0.1, 0.4, -0.2), q, m)
  expect_equal(r$mse, r$conditional_variance)
  s <- gaussian_target_mse(x, c(0.1, 0.4, -0.2), q, m, r$conditional_mean + 2)
  expect_equal(s$mse, r$mse + 4)
})
