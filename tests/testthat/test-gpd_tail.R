test_that("survival and quantile invert each other", {
  p <- c(0.9, 0.5, 0.1, 0.01, 1e-4)
  for (xi in c(-0.3, 0, 0.2, 0.7)) {
    y <- gpd_quantile_upper(p, scale = 2.5, shape = xi)
    expect_equal(gpd_survival(y, 2.5, xi), p, tolerance = 1e-9)
  }
})

test_that("density integrates to the survival difference", {
  for (xi in c(-0.3, 1e-4, 0.25)) {
    num <- stats::integrate(
      function(u) gpd_density(u, 2.5, xi),
      0.3,
      4
    )$value
    ana <- gpd_survival(0.3, 2.5, xi) - gpd_survival(4, 2.5, xi)
    expect_equal(num, ana, tolerance = 1e-7)
  }
})

test_that("a negative shape gives a finite upper endpoint", {
  xi <- -0.25
  scale <- 2
  endpoint <- -scale / xi
  expect_equal(gpd_survival(endpoint + 1, scale, xi), 0)
  expect_equal(gpd_density(endpoint + 1, scale, xi), 0)
  # Survival reaches zero continuously, so every positive probability is
  # still attainable strictly inside the support.
  expect_gt(gpd_survival(endpoint - 1e-6, scale, xi), 0)
})

test_that("the exponential limit is continuous at shape zero", {
  expect_equal(gpd_survival(3, 2, 1e-9), gpd_survival(3, 2, 0), tolerance = 1e-9)
  expect_equal(gpd_survival(3, 2, 0), exp(-1.5))
})

test_that("weighted fitting recovers known parameters", {
  set.seed(42)
  for (truth in list(c(1.5, 0.2), c(3, -0.15), c(2, 0.45))) {
    y <- gpd_quantile_upper(stats::runif(20000), truth[1], truth[2])
    fit <- fit_gpd_weighted(y)
    expect_equal(fit$scale, truth[1], tolerance = 0.15 * truth[1])
    expect_equal(fit$shape, truth[2], tolerance = 0.05)
  }
})

test_that("weights are equivalent to replicating points", {
  set.seed(7)
  y <- gpd_quantile_upper(stats::runif(300), 2, 0.2)
  doubled <- c(y, y[1:50])
  by_replication <- fit_gpd_weighted(doubled)
  by_weight <- fit_gpd_weighted(doubled, w = rep(1, 350))
  expect_equal(by_replication$shape, by_weight$shape, tolerance = 1e-4)
  expect_equal(by_replication$scale, by_weight$scale, tolerance = 1e-4)
})

test_that("shared-shape fitting recovers a common shape across varied scales", {
  set.seed(11)
  xi <- 0.18
  scales <- c(0.5, 1, 2, 4, 8, 16)
  excesses <- lapply(scales, function(s) {
    gpd_quantile_upper(stats::runif(3000), s, xi)
  })
  fit <- fit_gpd_shared_shape(excesses, n_boot = 20L)
  expect_equal(fit$shape, xi, tolerance = 0.03)
  expect_equal(fit$scale, scales, tolerance = 0.1)
  expect_true(fit$converged)
  # Plenty of information here, so the bootstrap spread should be tight.
  expect_lt(fit$shape_se, 0.02)
})

test_that("pooling beats the maximum of independent shape estimates", {
  # This is the failure mode the shared shape exists to fix: the equal-weight
  # mixture inherits max(xi_i), and that maximum is wildly inflated when each
  # xi_i comes from a small effective sample.
  set.seed(202)
  xi <- 0.18
  excesses <- lapply(stats::runif(100, 1, 5), function(s) {
    gpd_quantile_upper(stats::runif(15), s, xi)
  })
  independent <- vapply(
    excesses,
    function(y) fit_gpd_weighted(y)$shape,
    numeric(1L)
  )
  pooled <- fit_gpd_shared_shape(excesses, n_boot = 25L)

  expect_gt(max(independent), xi + 0.4)
  expect_lt(abs(pooled$shape - xi), 0.08)
  expect_lt(abs(pooled$shape - xi), abs(max(independent) - xi))
})

test_that("bootstrap bias correction moves the estimate toward the truth", {
  set.seed(99)
  xi <- 0.18
  make <- function() {
    lapply(stats::runif(100, 1, 5), function(s) {
      gpd_quantile_upper(stats::runif(15), s, xi)
    })
  }
  raw <- mean(replicate(6, fit_gpd_shared_shape(make(), n_boot = 0L)$shape))
  corrected <- mean(replicate(6, fit_gpd_shared_shape(
    make(),
    n_boot = 25L,
    bias_correct = TRUE
  )$shape))
  # The uncorrected pooled MLE is biased low (incidental parameters).
  expect_lt(raw, xi)
  expect_lt(abs(corrected - xi), abs(raw - xi))
})

test_that("degenerate input is handled without error", {
  expect_true(is.na(fit_gpd_shared_shape(list())$shape))
  expect_true(is.na(fit_gpd_shared_shape(list(1, numeric(0)))$shape))

  set.seed(5)
  good <- gpd_quantile_upper(stats::runif(50), 2, 0.2)
  fit <- fit_gpd_shared_shape(list(good, numeric(0), good), n_boot = 0L)
  expect_true(is.na(fit$scale[2]))
  expect_false(is.na(fit$scale[1]))

  # Non-positive and missing excesses are dropped rather than propagated.
  expect_true(is.finite(
    fit_gpd_shared_shape(list(c(good, -1, NA)), n_boot = 0L)$shape
  ))
})

test_that("the shape bias correction is off by default", {
  # Deliberate: correcting the shape alone makes return levels worse, because
  # scale noise inflates them in the opposite direction. See the function's
  # documentation and scripts/experiments/tail-index-pooling.R.
  set.seed(8)
  ex <- lapply(stats::runif(30, 1, 5), function(s) {
    gpd_quantile_upper(stats::runif(15), s, 0.15)
  })
  fit <- fit_gpd_shared_shape(ex, n_boot = 10L)
  expect_identical(fit$shape, fit$shape_raw)
  # ... and with n_boot > 0 the bootstrap still reports a standard error.
  expect_true(is.finite(fit$shape_se))
})

test_that("n_boot = 0 disables both the correction and the standard error", {
  set.seed(6)
  ex <- lapply(1:4, function(i) gpd_quantile_upper(stats::runif(100), i, 0.2))
  fit <- fit_gpd_shared_shape(ex, n_boot = 0L)
  expect_identical(fit$shape, fit$shape_raw)
  expect_true(is.na(fit$shape_se))
})
