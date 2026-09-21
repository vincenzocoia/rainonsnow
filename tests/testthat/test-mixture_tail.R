make_components <- function(seed = 3, k = 40) {
  set.seed(seed)
  list(
    threshold = stats::runif(k, 8, 12),
    tail_prob = stats::runif(k, 0.4, 0.6),
    scale = stats::runif(k, 1, 6),
    weights = rep(1 / k, k)
  )
}

# Deliberately naive reference implementation: one component at a time.
reference_survival <- function(cmp, shape, x) {
  vapply(
    x,
    function(xx) {
      sum(
        cmp$weights *
          cmp$tail_prob *
          gpd_survival(xx - cmp$threshold, cmp$scale, shape)
      )
    },
    numeric(1L)
  )
}

test_that("the vectorised survival matches a component-by-component sum", {
  cmp <- make_components()
  mt <- mixture_tail(cmp$threshold, cmp$tail_prob, cmp$scale, 0.17, cmp$weights)
  x <- seq(max(cmp$threshold), max(cmp$threshold) + 500, length.out = 50)
  expect_equal(
    mixture_tail_survival(mt, x),
    reference_survival(cmp, 0.17, x),
    tolerance = 1e-12
  )
  expect_true(all(diff(mixture_tail_survival(mt, x)) < 0))
})

test_that("evaluation below the graft region is refused, not extrapolated", {
  cmp <- make_components()
  mt <- mixture_tail(cmp$threshold, cmp$tail_prob, cmp$scale, 0.17, cmp$weights)
  expect_true(is.na(mixture_tail_survival(mt, mixture_tail_lower(mt) - 1)))
  expect_true(is.na(mixture_tail_quantile(mt, mixture_tail_max_prob(mt) * 1.01)))
})

test_that("quantile inverts survival", {
  cmp <- make_components()
  mt <- mixture_tail(cmp$threshold, cmp$tail_prob, cmp$scale, 0.17, cmp$weights)
  p <- c(0.1, 0.01, 1e-3, 1e-4, 1e-5)
  p <- p[p <= mixture_tail_max_prob(mt)]
  q <- mixture_tail_quantile(mt, p)
  expect_equal(mixture_tail_survival(mt, q), p, tolerance = 1e-9)
  expect_true(all(diff(q) > 0))
})

test_that("a single component reduces exactly to a shifted GP", {
  mt <- mixture_tail(10, 0.5, 3, 0.2)
  x <- c(10, 15, 40, 200)
  expect_equal(
    mixture_tail_survival(mt, x),
    0.5 * gpd_survival(x - 10, 3, 0.2),
    tolerance = 1e-14
  )
  # Repeating one component changes nothing.
  repeated <- mixture_tail(rep(10, 7), rep(0.5, 7), rep(3, 7), 0.2)
  expect_equal(
    mixture_tail_survival(repeated, x),
    mixture_tail_survival(mt, x),
    tolerance = 1e-14
  )
})

test_that("the single-GP approximation converges in the far tail", {
  cmp <- make_components()
  mt <- mixture_tail(cmp$threshold, cmp$tail_prob, cmp$scale, 0.17, cmp$weights)
  approx <- mixture_tail_gpd_approx(mt)
  x <- c(1e4, 1e7)
  rel <- abs(
    approx$tail_prob *
      gpd_survival(x - approx$threshold, approx$scale, approx$shape) -
      mixture_tail_survival(mt, x)
  ) /
    mixture_tail_survival(mt, x)
  expect_lt(rel[2], 1e-3)
  expect_lt(rel[2], rel[1])
})

test_that("a zero shape is handled", {
  cmp <- make_components()
  mt <- mixture_tail(cmp$threshold, cmp$tail_prob, cmp$scale, 0, cmp$weights)
  q <- mixture_tail_quantile(mt, 1e-4)
  expect_true(is.finite(q))
  expect_equal(mixture_tail_survival(mt, q), 1e-4, tolerance = 1e-10)
})

test_that("a bounded tail still reaches any positive probability", {
  cmp <- make_components()
  mt <- mixture_tail(cmp$threshold, cmp$tail_prob, cmp$scale, -0.3, cmp$weights)
  endpoint <- max(cmp$threshold + cmp$scale / 0.3)
  expect_equal(mixture_tail_survival(mt, endpoint + 1), 0)
  q <- mixture_tail_quantile(mt, 1e-12)
  expect_true(is.finite(q))
  expect_lt(q, endpoint)
})

test_that("unusable components are dropped and counted", {
  mt <- mixture_tail(
    threshold = c(10, NA, 12),
    tail_prob = c(0.5, 0.5, 0.5),
    scale = c(2, 3, -1),
    shape = 0.1
  )
  expect_length(mt$threshold, 1L)
  expect_equal(mt$n_dropped, 2L)
  expect_error(mixture_tail(NA_real_, 0.5, 1, 0.1))
})

test_that("scatter in the component shapes inflates the return level", {
  # The point of the shared shape. Independent per-hour fits have a spread of
  # roughly sd 0.3-0.4 at realistic tail sample sizes; that scatter alone
  # roughly doubles the 200-year level relative to a common shape.
  cmp <- make_components()
  shared <- mixture_tail_quantile(
    mixture_tail(cmp$threshold, cmp$tail_prob, cmp$scale, 0.15, cmp$weights),
    1 / 600
  )
  set.seed(1234)
  scattered <- median(replicate(200, {
    noisy <- pmin(pmax(stats::rnorm(40, 0.15, 0.3), -0.45), 1)
    mixture_tail_quantile(
      mixture_tail(cmp$threshold, cmp$tail_prob, cmp$scale, noisy, cmp$weights),
      1 / 600
    )
  }))
  expect_gt(scattered / shared, 1.3)
})

test_that("compression preserves the return curve", {
  set.seed(1)
  k <- 1000
  for (xi in c(0.05, 0.15, 0.35)) {
    mt <- mixture_tail(
      stats::runif(k, 8, 12),
      stats::runif(k, 0.4, 0.6),
      stats::runif(k, 1, 6),
      xi
    )
    small <- mixture_tail_compress(mt, 64L)
    expect_lte(length(small$threshold), 64L)
    # Far below the statistical uncertainty in any of these parameters.
    expect_lt(mixture_tail_compression_error(mt, small), 1e-3)
  }
})

test_that("compression is a no-op below the bin count and needs a shared shape", {
  mt <- mixture_tail(c(10, 11), c(0.5, 0.5), c(2, 3), 0.2)
  expect_identical(mixture_tail_compress(mt, 64L), mt)

  mixed <- mixture_tail(c(10, 11), c(0.5, 0.5), c(2, 3), c(0.1, 0.5))
  expect_error(mixture_tail_compress(mixed, 1L))
})

test_that("the density integrates to the survival difference", {
  cmp <- make_components()
  mt <- mixture_tail(cmp$threshold, cmp$tail_prob, cmp$scale, 0.17, cmp$weights)
  lo <- mixture_tail_lower(mt)
  num <- stats::integrate(
    function(u) mixture_tail_density(mt, u),
    lo,
    lo + 400,
    subdivisions = 500L
  )$value
  ana <- mixture_tail_survival(mt, lo) - mixture_tail_survival(mt, lo + 400)
  expect_equal(num, ana, tolerance = 1e-6)
  expect_true(is.na(mixture_tail_density(mt, lo - 1)))
})
