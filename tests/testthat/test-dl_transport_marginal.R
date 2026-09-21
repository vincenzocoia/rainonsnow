test_that("the copula parameter comes back from the ranks it was simulated with", {
  set.seed(20)
  rho <- 0.7
  n <- 4000L
  a <- stats::rnorm(n)
  b <- rho * a + sqrt(1 - rho^2) * stats::rnorm(n)
  expect_equal(dl_copula_par(a, b, "gaussian"), rho, tolerance = 0.05)
  # Clayton has tau = theta / (theta + 2), so a positive tau maps to a theta.
  th <- dl_copula_par(a, b, "survival_clayton")
  expect_equal(th / (th + 2), stats::cor(a, b, method = "kendall"),
               tolerance = 1e-8)
  # No Clayton has negative dependence.
  expect_true(is.na(dl_copula_par(a, -b, "clayton")))
  expect_true(is.na(dl_copula_par(1:3, 3:1, "gaussian")))
})

# A cell whose peak hours are the exact conditionals of a known process: hour i
# saw covariate x_i, so its predictive distribution is GP(scale = x_i, shape),
# grafted at zero. Transporting those through the process's own copula must
# return the process's marginal, from every hour.
exact_cell <- function(n = 200L, seed = 21) {
  set.seed(seed)
  d <- transport_dgp("pareto", alpha = 4, xi_cond = 0.1)
  s <- td_simulate(d, n)
  list(
    dgp = d,
    observed = s$y,
    df = data.frame(
      graft_of = 0,
      graft_tail_prob = 1,
      gp_scale = s$x,
      gp_shape = d$xi_cond
    )
  )
}

test_that("exact conditionals transport back to the exact marginal", {
  cell <- exact_cell()
  d <- cell$dgp
  me <- dl_cell_transport_ensemble(
    cell$df,
    family = function(sc, u) td_transport(d, sc, u),
    location = cell$df$gp_scale
  )
  expect_s3_class(me, "marginal_ensemble")
  expect_true(all(me$u > 0 & me$u < 1))

  # Every hour agrees, so the spread is one and the answer is the truth. The
  # ranks stand in for F_X(x), so this is right up to rank noise.
  expect_equal(dl_transport_spread(me, 1e-4), 1, tolerance = 0.25)
  expect_equal(
    marginal_ensemble_quantile(me, 1e-4) / td_marginal_quantile(d, 1e-4),
    1,
    tolerance = 0.1
  )
})

test_that("the wrong copula shows up as disagreement between hours", {
  cell <- exact_cell()
  d <- cell$dgp
  right <- dl_cell_transport_ensemble(
    cell$df,
    family = function(sc, u) td_transport(d, sc, u),
    location = cell$df$gp_scale
  )
  wrong <- dl_cell_transport_ensemble(
    cell$df,
    family = "gaussian",
    par = 0.2,
    location = cell$df$gp_scale
  )
  expect_lt(dl_transport_spread(right, 1e-4), 1.5)
  expect_gt(dl_transport_spread(wrong, 1e-4), 4)
})

test_that("the fitted route needs the observations, and drops unusable hours", {
  cell <- exact_cell()
  expect_error(dl_cell_transport_ensemble(cell$df), "`observed`")
  expect_error(
    dl_cell_transport_ensemble(cell$df, observed = cell$observed[1:5]),
    "one value per row"
  )
  expect_error(
    dl_cell_transport_ensemble(data.frame(graft_of = 1), observed = 1),
    "missing: gp_scale"
  )

  df <- cell$df
  df$gp_scale[1:20] <- NA
  me <- dl_cell_transport_ensemble(
    df,
    observed = cell$observed,
    family = "gaussian",
    location = cell$df$gp_scale
  )
  expect_equal(length(me$u), nrow(df) - 20L)

  # Too little to work with is NULL, not an error.
  expect_null(dl_cell_transport_ensemble(
    cell$df[1:3, ],
    observed = cell$observed[1:3],
    family = "gaussian",
    location = cell$df$gp_scale[1:3]
  ))
  # So is a predictive summary that does not vary: it carries no ranks, so
  # there is no copula to fit against it.
  expect_null(dl_cell_transport_ensemble(
    cell$df,
    observed = cell$observed,
    family = "gaussian",
    location = rep(1, nrow(cell$df))
  ))
})

test_that("return levels rise with the return period", {
  cell <- exact_cell()
  d <- cell$dgp
  me <- dl_cell_transport_ensemble(
    cell$df,
    family = function(sc, u) td_transport(d, sc, u),
    location = cell$df$gp_scale
  )
  rl <- dl_transport_return_levels(me, c(2, 10, 100, 1000), events_per_year = 4)
  expect_equal(nrow(rl), 4L)
  ok <- is.finite(rl$return_level)
  expect_true(all(diff(rl$return_level[ok]) > 0))
})
