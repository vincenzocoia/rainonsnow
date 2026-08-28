test_that("grid neighbours are the surrounding cells", {
  g <- expand.grid(x = 1:3, y = 1:3)
  nb <- grid_neighbours(g$x, g$y, radius = 1)
  # Centre of a 3x3 grid touches all eight others.
  centre <- which(g$x == 2 & g$y == 2)
  expect_length(nb[[centre]], 8L)
  expect_false(centre %in% nb[[centre]])
  # A corner touches three.
  corner <- which(g$x == 1 & g$y == 1)
  expect_length(nb[[corner]], 3L)
})

test_that("neighbours respect non-unit spacing and include_self", {
  g <- expand.grid(x = seq(0, 1, by = 0.25), y = seq(0, 1, by = 0.25))
  nb <- grid_neighbours(g$x, g$y, radius = 1)
  centre <- which(abs(g$x - 0.5) < 1e-9 & abs(g$y - 0.5) < 1e-9)
  expect_length(nb[[centre]], 8L)

  with_self <- grid_neighbours(g$x, g$y, radius = 1, include_self = TRUE)
  expect_length(with_self[[centre]], 9L)
  expect_true(centre %in% with_self[[centre]])
})

test_that("smoothing reduces error when the field is smooth", {
  set.seed(1)
  g <- expand.grid(x = 1:6, y = 1:6)
  truth <- 0.1 + 0.02 * g$x
  se <- rep(0.15, nrow(g))

  raw_err <- smooth_err <- numeric(30)
  for (i in seq_along(raw_err)) {
    est <- truth + stats::rnorm(nrow(g), sd = se)
    out <- smooth_tail_shape(est, se, g$x, g$y)
    raw_err[i] <- mean((est - truth)^2)
    smooth_err[i] <- mean((out$shape - truth)^2)
  }
  expect_lt(mean(smooth_err), mean(raw_err))
})

test_that("shrinkage tracks each cell's own precision", {
  # A real shape field: enough genuine between-cell variation that tau^2 > 0,
  # with precision differing across cells.
  set.seed(3)
  g <- expand.grid(x = 1:6, y = 1:6)
  # A trend PLUS genuine cell-to-cell roughness, so there is real between-cell
  # variation for the smoother to preserve (tau^2 > 0).
  truth <- 0.05 + 0.05 * g$x + stats::rnorm(nrow(g), sd = 0.08)
  se <- rep(0.05, nrow(g))
  se[c(10, 20)] <- 0.5 # two badly determined cells
  est <- truth + stats::rnorm(nrow(g), sd = se)

  out <- smooth_tail_shape(est, se, g$x, g$y)
  expect_gt(out$tau, 0)

  # Precise cells keep much more of their own estimate than imprecise ones.
  # The absolute weight depends on how tau^2 comes out; the ordering is the
  # property that has to hold.
  expect_gt(mean(out$weight[-c(10, 20)]), 5 * max(out$weight[c(10, 20)]))
  expect_lt(max(out$weight[c(10, 20)]), 0.2)

  # And the badly determined cells end up much closer to the truth.
  expect_lt(
    mean(abs(out$shape[c(10, 20)] - truth[c(10, 20)])),
    mean(abs(est[c(10, 20)] - truth[c(10, 20)]))
  )
})

test_that("neighbours are combined by precision, not equally", {
  # One neighbour is far off but barely determined; it should not dominate.
  g <- expand.grid(x = 1:3, y = 1:3)
  est <- rep(0.15, 9)
  est[1] <- 1.0
  se <- rep(0.05, 9)
  se[1] <- 10 # that outlier carries almost no information

  centre <- which(g$x == 2 & g$y == 2)
  out <- smooth_tail_shape(est, se, g$x, g$y)
  expect_lt(abs(out$shape[centre] - 0.15), 0.05)
})

test_that("indistinguishable shapes collapse to a common value", {
  # All the spread is sampling noise, so tau^2 estimates to zero and every cell
  # should take its neighbourhood value rather than its own.
  set.seed(2)
  g <- expand.grid(x = 1:5, y = 1:5)
  est <- 0.15 + stats::rnorm(nrow(g), sd = 0.05)
  out <- smooth_tail_shape(est, se = rep(1, nrow(g)), x = g$x, y = g$y)
  expect_equal(out$tau, 0)
  expect_true(all(out$weight == 0))
  expect_lt(stats::sd(out$shape), stats::sd(est))
})

test_that("missing and degenerate input is passed through", {
  g <- expand.grid(x = 1:3, y = 1:3)
  est <- c(rep(0.1, 8), NA)
  out <- smooth_tail_shape(est, se = rep(0.1, 9), x = g$x, y = g$y)
  expect_true(is.na(out$shape[9]))

  one <- smooth_tail_shape(0.2, 0.1, 1, 1)
  expect_equal(one$shape, 0.2)
})
