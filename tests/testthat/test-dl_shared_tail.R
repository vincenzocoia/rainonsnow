# These need distionary, which supplies the conditional quantiles and survival
# probabilities the standardised fit reads off each predictive distribution.
skip_if_not_installed("distionary")

make_cell <- function(k = 400L, xi = 0.2, seed = 5) {
  set.seed(seed)
  scale <- stats::runif(k, 1, 5)
  base <- stats::runif(k, 5, 9)
  observed <- base + gpd_quantile_upper(stats::runif(k), scale, xi)
  support <- sort(c(base, base + gpd_quantile_upper(stats::runif(k), scale, xi)))
  distributions <- lapply(seq_len(k), function(i) {
    w <- stats::dnorm(support, mean = base[i] + scale[i], sd = 2.5 * scale[i])
    distionary::dst_empirical(support, weights = w / sum(w))
  })
  list(distributions = distributions, observed = observed, xi = xi)
}

test_that("the standardised fit recovers a shared shape", {
  cell <- make_cell()
  fit <- fit_shared_tail_standardised(cell$distributions, cell$observed)
  expect_equal(fit$gp_shape, cell$xi, tolerance = 0.25)
  expect_gt(fit$n_exceedances, 50L)
  expect_true(all(is.finite(fit$graft_of)))
  expect_true(all(fit$gp_scale > 0, na.rm = TRUE))
})

test_that("the back-transform is a pure rescaling", {
  # y = loc + sc * z, so an excess in z maps to sc * excess in y: the shape is
  # unchanged and every hour's scale is its own gap times one common factor.
  cell <- make_cell()
  # State the quantiles here rather than relying on the defaults, so the test
  # checks the back-transform and not the current default choice.
  fit <- fit_shared_tail_standardised(
    cell$distributions,
    cell$observed,
    p_location = 0.5,
    p_scale = 0.9
  )
  gap <- vapply(
    cell$distributions,
    function(d) {
      distionary::eval_quantile(d, at = 0.9) - distionary::eval_quantile(d, at = 0.5)
    },
    numeric(1L)
  )
  ratio <- fit$gp_scale / gap
  expect_equal(stats::sd(ratio), 0, tolerance = 1e-9)
})

test_that("each hour gets its own graft point", {
  cell <- make_cell()
  fit <- fit_shared_tail_standardised(cell$distributions, cell$observed)
  # A common threshold on the standardised scale is a different value in each
  # hour's own units.
  expect_gt(stats::sd(fit$graft_of, na.rm = TRUE), 0)
})

test_that("the standard error is the usual asymptotic one", {
  cell <- make_cell()
  fit <- fit_shared_tail_standardised(cell$distributions, cell$observed)
  # Available because the pooled values are one per observation rather than a
  # re-weighting of the same responses many times over.
  expect_equal(
    fit$gp_shape_se,
    sqrt((1 + fit$gp_shape)^2 / fit$n_exceedances),
    tolerance = 1e-12
  )
})

test_that("the wrapper routes to the standardised fit and stays compatible", {
  cell <- make_cell()
  direct <- fit_shared_tail_standardised(cell$distributions, cell$observed)
  via <- dl_fit_cell_shared_tail(
    cell$distributions,
    method = "standardise",
    observed = cell$observed
  )
  expect_equal(via$gp_shape, direct$gp_shape)
  expect_true(all(via$graft_tail_prob > 0 & via$graft_tail_prob < 1, na.rm = TRUE))

  # The output still feeds the closed-form mixture tail unchanged.
  df <- data.frame(
    graft_of = via$graft_of,
    graft_tail_prob = via$graft_tail_prob,
    gp_scale = via$gp_scale,
    gp_shape = via$gp_shape
  )
  expect_true(dl_cell_mixture_tail(df)$shared_shape)
})

test_that("the standardised route needs the observations", {
  cell <- make_cell()
  expect_error(
    dl_fit_cell_shared_tail(cell$distributions, method = "standardise"),
    "observed"
  )
  expect_error(
    fit_shared_tail_standardised(cell$distributions, cell$observed[1:10]),
    "one value per distribution"
  )
})

test_that("a supplied shape still uses the pooled branch", {
  # Passing `shape` means the caller already has a shape and wants the scales
  # refitted under it -- as the animation app does with a cell's own shape.
  cell <- make_cell()
  fit <- dl_fit_cell_shared_tail(
    cell$distributions,
    method = "standardise",
    observed = cell$observed,
    shape = 0.3
  )
  expect_equal(fit$gp_shape, 0.3)
})

test_that("too little usable data returns the empty result", {
  cell <- make_cell()
  fit <- fit_shared_tail_standardised(
    cell$distributions[1:5],
    cell$observed[1:5]
  )
  expect_true(is.na(fit$gp_shape))
  expect_equal(fit$n_exceedances, 0L)
})
