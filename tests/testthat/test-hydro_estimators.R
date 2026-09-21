test_that("the GEV fit recovers what it was simulated from", {
  set.seed(3)
  n <- 4000L
  y <- 10 + 2 / 0.15 * ((-base::log(stats::runif(n)))^-0.15 - 1)
  f <- fit_gev(y)
  expect_equal(f$location, 10, tolerance = 0.1)
  expect_equal(f$scale, 2, tolerance = 0.15)
  expect_equal(f$shape, 0.15, tolerance = 0.05)
  # The 100-year level, against the closed form.
  expect_equal(
    gev_return_level(f, 0.01),
    10 + 2 / 0.15 * ((-base::log(0.99))^-0.15 - 1),
    tolerance = 0.05
  )
})

test_that("the generalized Pareto fit recovers its parameters", {
  set.seed(4)
  x <- gpd_quantile_upper(stats::runif(5000L), 2, 0.2)
  f <- fit_gp_ml(x)
  expect_equal(f$scale, 2, tolerance = 0.15)
  expect_equal(f$shape, 0.2, tolerance = 0.06)
})

test_that("the structural fit recovers the two drivers from the max alone", {
  # Neither driver is observed: only the covariate and the peak.
  set.seed(5)
  d <- two_process_dgp()
  s <- tp_simulate(d, 5000L)
  f <- fit_two_process(s$x, s$y)
  expect_equal(f$scale_a, d$scale_a, tolerance = 0.12)
  expect_equal(f$shape_a, d$shape_a, tolerance = 0.08)
  expect_equal(f$b_coef, d$b_coef, tolerance = 0.12)
  expect_equal(f$shape_b, d$shape_b, tolerance = 0.08)
  expect_true(f$shape_a <= f$shape_b)
})

test_that("the scaled fit beats a marginal fit for a covariate-scaled driver", {
  # B's marginal is a scale mixture, so a generalized Pareto fitted to it is
  # misspecified and reads too heavy; fitted on the covariate it is exact.
  set.seed(6)
  d <- two_process_dgp()
  s <- tp_simulate(d, 5000L)
  cond <- fit_scaled_gp(s$x, s$b)
  expect_equal(cond$b_coef, d$b_coef, tolerance = 0.12)
  expect_equal(cond$shape, d$shape_b, tolerance = 0.06)
  expect_gt(fit_gp_ml(s$b)$shape, d$shape_b + 0.05)
})

test_that("averaging over the observed covariate collapses far out, but not near", {
  # Both use the same conditional model, so the only difference is which
  # covariate distribution the conditional is removed over. Averaging over the
  # OBSERVED values truncates that distribution at the sample maximum, and once
  # the bounded driver has died the average has nothing left to give.
  #
  # Where that bites is the point. At the levels a hydrologist asks about --
  # the 20- to 200-year event, y of 3 to 11 here -- the two agree to within
  # about ten percent. The truncation only takes over far beyond them.
  set.seed(7)
  d <- two_process_dgp()
  s <- tp_simulate(d, 400L)
  f <- fit_two_process(s$x, s$y)
  fx <- fit_gp_ml(s$x)

  design <- tp_marginal_quantile(d, 1 / c(20, 200))
  near <- tpf_marginal_survival(f, fx, design) /
    tpf_mixture_survival(f, s$x, design)
  expect_true(all(near > 0.8 & near < 1.25))

  far <- 500
  expect_gt(
    tpf_marginal_survival(f, fx, far) / tpf_mixture_survival(f, s$x, far),
    100
  )
})

test_that("the structural marginal agrees with the truth given true parameters", {
  d <- two_process_dgp()
  f <- list(scale_a = d$scale_a, shape_a = d$shape_a,
            b_coef = d$b_coef, shape_b = d$shape_b)
  fx <- list(scale = d$scale_x, shape = d$shape_x)
  y <- c(3, 8, 25)
  expect_equal(tpf_marginal_survival(f, fx, y), tp_marginal_survival(d, y),
               tolerance = 1e-3)
})

test_that("invert_survival inverts, and declines when it cannot", {
  expect_equal(invert_survival(function(y) exp(-y), c(0.5, 0.01)),
               c(base::log(2), base::log(100)), tolerance = 1e-6)
  # A survival that never gets below the target has no answer to give.
  expect_true(is.na(invert_survival(function(y) 0.5, 0.01)))
})

test_that("the upper-region copula fit recovers the parameter it was given", {
  set.seed(8)
  n <- 3000L
  u <- stats::runif(n)
  v <- cop_h_inverse(stats::runif(n), u, "gaussian", 0.6)
  expect_equal(fit_copula_upper(u, v, "gaussian", 0.5), 0.6, tolerance = 0.08)
  expect_true(is.na(fit_copula_upper(u[1:4], v[1:4], "gaussian")))
})

test_that("the copula density integrates to one", {
  g <- seq(0.004, 0.996, length.out = 250L)
  for (fam in c("gaussian", "clayton", "survival_clayton")) {
    par <- if (fam == "gaussian") 0.6 else 2
    m <- outer(g, g, function(a, b) cop_density(a, b, fam, par))
    expect_equal(mean(m), 1, tolerance = 0.02, info = fam)
  }
})
