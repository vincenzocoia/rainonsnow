test_that("the bounded driver really is bounded, and the peak is their max", {
  d <- two_process_dgp(scale_a = 1, shape_a = -0.25)
  expect_equal(d$endpoint_a, 4)
  expect_equal(gpd_survival(4.001, d$scale_a, d$shape_a), 0)
  set.seed(1)
  s <- tp_simulate(d, 500L)
  expect_equal(s$y, pmax(s$a, s$b))
  expect_true(max(s$a) < d$endpoint_a)
})

test_that("the conditional survival is the competing-risk one", {
  d <- two_process_dgp()
  y <- c(1, 3, 8)
  x <- 2
  sa <- gpd_survival(y, d$scale_a, d$shape_a)
  sb <- gpd_survival(y, d$b_coef * x, d$shape_b)
  expect_equal(tp_conditional_survival(d, y, x), 1 - (1 - sa) * (1 - sb))
})

test_that("the cached marginal matches the quadrature and simulation", {
  d <- two_process_dgp()
  y <- c(1, 3, 8, 30, 200)
  expect_equal(
    tp_marginal_survival(d, y),
    tp_marginal_survival(d, y, ngrid = 30001L, exact = TRUE),
    tolerance = 1e-3
  )
  set.seed(2)
  s <- tp_simulate(d, 2e5L)
  for (yy in c(1, 3, 8)) {
    expect_equal(mean(s$y > yy), tp_marginal_survival(d, yy), tolerance = 0.05)
  }
})

test_that("the quantile inverts the survival", {
  d <- two_process_dgp()
  p <- 1 / c(2, 20, 200, 2000)
  expect_equal(tp_marginal_survival(d, tp_marginal_quantile(d, p)), p,
               tolerance = 1e-3)
})

test_that("the marginal is kinked, which is the point of the construction", {
  # Below the bounded driver's endpoint the light component still contributes;
  # above it, only the heavy one is left. So the local index is far from its
  # asymptote across the whole design range, and a single generalized Pareto
  # cannot describe both halves.
  d <- two_process_dgp()
  q <- tp_marginal_quantile(d, 1 / c(20, 200))
  evi <- tp_local_evi(d, q)
  expect_true(all(evi > 0.4))
  expect_gt(max(evi) / max(d$shape_x, d$shape_b), 1.7)
})

test_that("the covariate selects the regime, not just the scale", {
  # This is what a generic copula cannot do: at low rainfall the conditional is
  # the bounded driver's, at high rainfall the heavy driver's.
  d <- two_process_dgp()
  lo <- tp_x_quantile(d, 0.2)
  hi <- tp_x_quantile(d, 0.95)
  y <- c(3.5, 3.9)
  # Near the bounded driver's endpoint the low-rainfall conditional has all but
  # died while the high-rainfall one has not.
  expect_lt(
    tp_conditional_survival(d, y[2], lo) / tp_conditional_survival(d, y[1], lo),
    tp_conditional_survival(d, y[2], hi) / tp_conditional_survival(d, y[1], hi)
  )
})

test_that("the true transport returns the marginal from every covariate value", {
  d <- two_process_dgp()
  u <- c(0.1, 0.5, 0.9)
  y0 <- 8
  s_cond <- tp_conditional_survival(d, y0, tp_x_quantile(d, u))
  expect_equal(
    tp_transport(d, s_cond, u),
    rep(tp_marginal_survival(d, y0), 3L),
    tolerance = 1e-4
  )
})

test_that("the DGP insists on a bounded driver", {
  expect_error(two_process_dgp(shape_a = 0.1), "must be negative")
  expect_error(two_process_dgp(scale_a = -1), "positive")
})
