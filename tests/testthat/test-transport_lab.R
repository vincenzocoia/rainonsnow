test_that("the derived marginal is a survival function", {
  d <- transport_dgp()
  expect_equal(td_marginal_survival(d, 0), 1)
  y <- exp(seq(-3, 20, length.out = 60L))
  s <- td_marginal_survival(d, y)
  expect_true(all(s > 0 & s <= 1))
  expect_true(all(diff(s) < 0))
})

test_that("the cached curve matches the quadrature it was built from", {
  for (d in list(
    transport_dgp(),
    transport_dgp(sdlog = 2),
    transport_dgp("pareto", alpha = 4, xi_cond = 0.1),
    transport_dgp("pareto", alpha = 8)
  )) {
    y <- 10^c(0, 2, 4, 8, 16)
    expect_equal(
      td_marginal_survival(d, y),
      td_marginal_survival(d, y, ngrid = 60001L, exact = TRUE),
      tolerance = 1e-3,
      info = d$x_dist
    )
  }
})

test_that("the quadrature agrees with simulation where simulation can reach", {
  set.seed(9)
  d <- transport_dgp("pareto", alpha = 4, xi_cond = 0.1)
  s <- td_simulate(d, 3e5L)
  for (y in c(1, 5, 20)) {
    expect_equal(mean(s$y > y), td_marginal_survival(d, y), tolerance = 0.05)
  }
})

test_that("the quantile inverts the survival", {
  d <- transport_dgp("pareto", alpha = 4)
  p <- c(1e-2, 1e-5, 1e-10)
  expect_equal(td_marginal_survival(d, td_marginal_quantile(d, p)), p,
               tolerance = 1e-3)
})

test_that("the local tail index reaches the asymptote Breiman gives it", {
  # Pareto covariate heavier than the conditional: the marginal takes the
  # covariate's index, and CEVI is the ratio.
  d <- transport_dgp("pareto", alpha = 4, xi_cond = 0.1)
  expect_equal(td_asymptotic_evi(d), 0.25)
  expect_equal(td_local_evi(d, td_marginal_quantile(d, 1e-9)), 0.25,
               tolerance = 1e-3)
  # ... but it is still far from it at anything a sample reaches, which is what
  # makes a canned generalized Pareto fit biased rather than merely imprecise.
  expect_gt(td_local_evi(d, td_marginal_quantile(d, 1e-3)), 0.26)

  # Lognormal covariate: every moment is finite, so the marginal inherits the
  # CONDITIONAL index and CEVI is one.
  d2 <- transport_dgp("lognormal", sdlog = 1, xi_cond = 0.2)
  expect_equal(td_asymptotic_evi(d2), 0.2)
  expect_gt(td_local_evi(d2, td_marginal_quantile(d2, 1e-3)), 0.35)
})

test_that("the conditional is exactly generalized Pareto", {
  d <- transport_dgp()
  y <- c(1, 10, 100)
  expect_equal(
    td_conditional_survival(d, y, x = 3),
    gpd_survival(y, scale = 3, shape = d$xi_cond)
  )
})

test_that("the true transport returns the same marginal from every covariate value", {
  # This is the identity the whole method rests on: with the true copula and the
  # true conditional, every x transports to the same answer.
  d <- transport_dgp("pareto", alpha = 4, xi_cond = 0.1)
  u <- c(0.05, 0.3, 0.6, 0.95)
  y0 <- 300
  s_cond <- td_conditional_survival(d, y0, td_x_quantile(d, u))
  expect_equal(
    td_transport(d, s_cond, u),
    rep(td_marginal_survival(d, y0), length(u)),
    tolerance = 1e-6
  )
})

test_that("h and the transport are inverses", {
  d <- transport_dgp()
  s <- c(1e-2, 1e-5, 1e-9)
  expect_equal(td_transport(d, td_h_surv(d, s, 0.8), 0.8), s, tolerance = 1e-3)
})

test_that("the covariate distribution is the one asked for", {
  d <- transport_dgp("pareto", alpha = 4)
  expect_equal(td_x_cdf(d, td_x_quantile(d, c(0.1, 0.9))), c(0.1, 0.9))
  expect_equal(td_x_quantile(d, 0), 1)
  d2 <- transport_dgp("lognormal", sdlog = 1.5)
  expect_equal(td_x_cdf(d2, td_x_quantile(d2, 0.37)), 0.37)
})

test_that("the DGP rejects impossible parameters", {
  expect_error(transport_dgp(xi_cond = 0), "xi_cond")
  expect_error(transport_dgp(sdlog = -1), "positive")
})

test_that("learned conditionals have the shape the tail machinery expects", {
  set.seed(12)
  d <- transport_dgp("pareto", alpha = 4, xi_cond = 0.1)
  s <- td_simulate(d, 2000L)
  cc <- td_conditionals(s$u, s$y, n_anchor = 12L, bandwidth = 0.05)
  expect_equal(length(cc$pieces), 12L)
  expect_equal(length(cc$u), 12L)
  p <- cc$pieces[[6]]
  expect_true(all(c("threshold", "excess", "weight", "tail_prob", "n_eff") %in%
                    names(p)))
  expect_true(all(p$excess > 0))
  expect_equal(p$tail_prob, sum(p$weight))
  # A conditional at a high covariate value sits above one at a low value.
  expect_gt(cc$pieces[[12]]$threshold, cc$pieces[[1]]$threshold)
})

test_that("the learned conditionals recover the conditional shape", {
  set.seed(13)
  d <- transport_dgp("pareto", alpha = 4, xi_cond = 0.1)
  s <- td_simulate(d, 5000L)
  cc <- td_conditionals(s$u, s$y, n_anchor = 30L, bandwidth = 0.03)
  me <- transport_ensemble(
    cc$pieces,
    cc$u,
    family = function(sc, u) td_transport(d, sc, u)
  )
  expect_equal(me$shape[1], d$xi_cond, tolerance = 0.05)
  # And the combined estimate lands near the derived truth, which a generalized
  # Pareto fitted to Y alone does not.
  truth <- td_marginal_quantile(d, 1e-5)
  expect_equal(marginal_ensemble_quantile(me, 1e-5) / truth, 1,
               tolerance = 0.15)
})

test_that("the canned fit reads the local slope, not the asymptotic one", {
  set.seed(14)
  d <- transport_dgp("lognormal", sdlog = 1, xi_cond = 0.2)
  g <- td_canned_gp(td_simulate(d, 5000L)$y)
  expect_gt(g$shape, 0.3)
  expect_equal(g$tail_prob, 0.1, tolerance = 1e-3)
  expect_true(all(diff(g$survival(c(20, 200, 2000))) < 0))
})
