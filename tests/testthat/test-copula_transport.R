test_that("h and its inverse round-trip in probability space", {
  for (fam in c("gaussian", "clayton", "survival_clayton")) {
    par <- if (fam == "gaussian") 0.7 else 1.5
    v <- c(0.01, 0.3, 0.7, 0.99)
    u <- c(0.05, 0.4, 0.6, 0.95)
    expect_equal(cop_h_inverse(cop_h(v, u, fam, par), u, fam, par), v,
                 tolerance = 1e-8, info = fam)
  }
})

test_that("survival space agrees with probability space where both work", {
  s <- c(1e-1, 1e-2, 1e-4)
  for (fam in c("gaussian", "clayton", "survival_clayton")) {
    par <- if (fam == "gaussian") 0.7 else 1.5
    expect_equal(
      cop_h_surv(s, 0.6, fam, par),
      1 - cop_h(1 - s, 0.6, fam, par),
      tolerance = 1e-6,
      info = fam
    )
  }
})

test_that("survival space survives the far tail, where the naive route does not", {
  # 1 - h^{-1}(1 - s) forms 1 - s for a tiny s, which rounds to one. At 1e-10
  # the round trip through probability space comes back about 130 times too
  # large; through survival space it is exact.
  s <- 1e-10
  for (fam in c("gaussian", "clayton", "survival_clayton")) {
    par <- if (fam == "gaussian") 0.7 else 1.5
    back <- cop_h_surv_inverse(cop_h_surv(s, 0.4, fam, par), 0.4, fam, par)
    expect_equal(back, s, tolerance = 1e-9, info = fam)
  }
  naive <- 1 - cop_h_inverse(1 - cop_h_surv(s, 0.4, "gaussian", 0.7), 0.4,
                             "gaussian", 0.7)
  expect_gt(naive / s, 10)
})

test_that("the Gaussian h-function is the textbook one", {
  v <- 0.8
  u <- 0.3
  rho <- 0.6
  expect_equal(
    cop_h(v, u, "gaussian", rho),
    stats::pnorm((stats::qnorm(v) - rho * stats::qnorm(u)) / sqrt(1 - rho^2))
  )
})

test_that("the Clayton h-function matches its direct formula", {
  v <- 0.4
  u <- 0.7
  th <- 2
  expect_equal(
    cop_h(v, u, "clayton", th),
    u^(-th - 1) * (u^-th + v^-th - 1)^(-1 / th - 1)
  )
})

test_that("a function may stand in for a copula family", {
  f <- function(s, u) s * u
  expect_equal(cop_h_surv_inverse(0.01, 0.5, f), 0.005)
  expect_equal(transport_survival(0.01, 0.5, f), 0.005)
})

make_ensemble <- function(...) {
  marginal_ensemble(
    u = c(0.2, 0.5, 0.8),
    threshold = c(9, 10, 11),
    tail_prob = rep(0.3, 3),
    scale = c(2, 2.5, 3),
    shape = 0.1,
    family = "gaussian",
    par = 0.7,
    ...
  )
}

test_that("an ensemble drops unusable conditionals and keeps the rest", {
  me <- marginal_ensemble(
    u = c(0.2, 1.5, 0.8),
    threshold = c(9, 10, 11),
    tail_prob = c(0.3, 0.3, NA),
    scale = 2,
    shape = 0.1
  )
  expect_equal(length(me$u), 1L)
  expect_equal(me$n_dropped, 2L)
  expect_error(
    marginal_ensemble(0.5, 10, 0.3, -1, 0.1),
    "No usable conditionals"
  )
})

test_that("every conditional gives its own estimate of the same marginal", {
  me <- make_ensemble()
  paths <- marginal_ensemble_paths(me, c(20, 50, 200))
  expect_equal(dim(paths), c(3L, 3L))
  # Estimates of one quantity, so they are the same order of magnitude but not
  # identical -- this ensemble's conditionals disagree by construction.
  expect_true(all(paths > 0 & paths < 1))
  expect_true(all(diff(paths[, 1]) < 0))
})

test_that("the combined curve is a survival function", {
  me <- make_ensemble()
  y <- exp(seq(log(12), log(5000), length.out = 40L))
  for (cm in c("median", "mean", "weighted_mean")) {
    s <- marginal_ensemble_survival(me, y, cm)
    expect_true(all(s > 0 & s < 1), info = cm)
    expect_true(all(diff(s) < 0), info = cm)
  }
})

test_that("the median is not dragged by one heavy conditional", {
  # A single rogue shape moves the mean a long way and the median by one rank.
  base <- marginal_ensemble(
    u = seq(0.1, 0.9, length.out = 9L),
    threshold = 10,
    tail_prob = 0.3,
    scale = 2,
    shape = 0.1,
    par = 0.7
  )
  rogue <- base
  rogue$shape[5] <- 0.9
  y <- 500
  expect_equal(
    marginal_ensemble_survival(rogue, y, "median"),
    marginal_ensemble_survival(base, y, "median"),
    tolerance = 0.15
  )
  expect_gt(
    marginal_ensemble_survival(rogue, y, "mean") /
      marginal_ensemble_survival(base, y, "mean"),
    5
  )
})

test_that("the quantile inverts the combined survival", {
  me <- make_ensemble()
  p <- c(1e-2, 1e-4, 1e-6)
  for (cm in c("median", "mean", "weighted_mean")) {
    q <- marginal_ensemble_quantile(me, p, cm)
    expect_equal(marginal_ensemble_survival(me, q, cm), p,
                 tolerance = 1e-6, info = cm)
  }
  # Above what the ensemble can describe there is no answer to give.
  expect_true(is.na(marginal_ensemble_quantile(me, 0.99)))
})

test_that("transport_ensemble shares one shape unless told not to", {
  set.seed(4)
  pieces <- lapply(1:6, function(i) {
    list(
      threshold = 10 * i,
      excess = stats::rexp(80L, 1 / i),
      weight = rep(1 / 80, 80L),
      tail_prob = 0.3,
      n_eff = 80
    )
  })
  u <- seq(0.15, 0.85, length.out = 6L)

  shared <- transport_ensemble(pieces, u, par = 0.7)
  expect_equal(length(unique(shared$shape)), 1L)
  expect_equal(length(unique(shared$scale)), 6L)

  free <- transport_ensemble(pieces, u, par = 0.7, shared_shape = FALSE)
  expect_gt(stats::sd(free$shape), 0)

  fixed <- transport_ensemble(pieces, u, par = 0.7, shape = 0.25)
  expect_equal(unique(fixed$shape), 0.25)

  # Weights default to each conditional's effective sample size.
  expect_equal(shared$weights, rep(1 / 6, 6L))

  # NULL entries are dropped, and the u values follow them.
  pieces[2] <- list(NULL)
  dropped <- transport_ensemble(pieces, u, par = 0.7)
  expect_equal(length(dropped$u), 5L)
  expect_false(0.29 %in% round(dropped$u, 2))

  expect_error(transport_ensemble(pieces, u[1:3]), "same length")
})
