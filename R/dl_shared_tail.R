# Bridge between the peak-hour predictive distributions produced by script 4 and
# the dependency-free tail machinery in gpd_tail.R / mixture_tail.R.
#
# Only `dl_tail_pieces()` touches distionary. Everything downstream of it works
# on plain numeric vectors, which is what keeps the marginal and return-level
# steps fast and testable.

#' Extract the upper tail of an empirical predictive distribution
#'
#' Splits a finite (empirical) distribution at the same adaptive threshold
#' [fit_and_graft_gp()] uses: the smaller of the weighted quantile and the
#' unweighted quantile of the support, so a tail carried by a handful of
#' heavily-weighted outcomes does not drag the threshold up.
#'
#' @param dst A finite distribution, as returned by `distionary::dst_empirical()`
#'   (the output of [predict.dl_rqforest()]).
#' @param adaptive_threshold Number between 0 and 1; see [fit_and_graft_gp()].
#' @returns A list with `threshold`, `excess`, `weight`, `tail_prob` (the mass
#'   above the threshold) and `n_eff` (the Kish effective sample size of the
#'   exceedances). Returns `NULL` for a null or unusable distribution.
#' @examples
#' \dontrun{
#' d <- predict(model)[[1]]
#' dl_tail_pieces(d)
#' }
#' @export
dl_tail_pieces <- function(dst, adaptive_threshold = 0.5) {
  if (is.null(dst) || !inherits(dst, "dst")) {
    return(NULL)
  }
  pars <- try(distionary::parameters(dst), silent = TRUE)
  if (inherits(pars, "try-error")) {
    return(NULL)
  }
  outcomes <- pars$outcomes
  probs <- pars$probs
  if (is.null(outcomes) || is.null(probs) || length(outcomes) < 3L) {
    return(NULL)
  }

  # Same threshold rule as fit_and_graft_gp(), computed directly from the
  # support so this needs neither distplyr::trim_left() nor a second pass.
  weighted_q <- try(
    distionary::eval_quantile(dst, at = adaptive_threshold),
    silent = TRUE
  )
  unweighted_q <- stats::quantile(
    outcomes,
    adaptive_threshold,
    names = FALSE,
    na.rm = TRUE
  )
  threshold <- if (inherits(weighted_q, "try-error") || !is.finite(weighted_q)) {
    unweighted_q
  } else {
    min(weighted_q, unweighted_q)
  }

  above <- is.finite(outcomes) & is.finite(probs) & outcomes > threshold
  if (sum(above) < 2L) {
    return(NULL)
  }
  w <- probs[above]
  list(
    threshold = threshold,
    excess = outcomes[above] - threshold,
    weight = w,
    tail_prob = sum(w),
    n_eff = sum(w)^2 / sum(w^2)
  )
}

#' Fit one shared GP tail shape across a cell's peak hours
#'
#' Replaces the per-hour independent fits of [fit_and_graft_gp()] with a single
#' shape for the whole cell and a scale per hour. This is what stops the cell's
#' marginal inheriting `max_i(xi_i)`; see the note at the top of
#' `R/mixture_tail.R` and the experiment in
#' `scripts/experiments/tail-index-pooling.R`.
#'
#' @param distributions A list of predictive distributions for one cell (the
#'   `distribution_forest` column for that cell).
#' @param adaptive_threshold Passed to [dl_tail_pieces()].
#' @param n_boot Bootstrap replicates for the shape standard error (and the
#'   bias correction when enabled); see [fit_gpd_shared_shape()].
#' @param bias_correct Whether to bias-correct the shape. Off by default --
#'   see the \dQuote{Why the correction is off by default} section of
#'   [fit_gpd_shared_shape()].
#' @param shape Optionally fix the shape instead of estimating it -- used to
#'   apply a spatially smoothed shape from [smooth_tail_shape()] while still
#'   refitting each hour's scale.
#' @returns A list with per-hour vectors `graft_of`, `graft_tail_prob`,
#'   `gp_scale`, and the scalars `gp_shape`, `gp_shape_se`, `n_hours_used`,
#'   `n_eff_median`.
#' @seealso [dl_encode_peak_hour_distributions()], [mixture_tail()]
#' @export
dl_fit_cell_shared_tail <- function(
  distributions,
  adaptive_threshold = 0.5,
  n_boot = 25L,
  shape = NULL,
  bias_correct = FALSE
) {
  pieces <- lapply(distributions, dl_tail_pieces, adaptive_threshold)
  usable <- !vapply(pieces, is.null, logical(1L))
  n <- length(distributions)

  empty <- list(
    graft_of = rep(NA_real_, n),
    graft_tail_prob = rep(NA_real_, n),
    gp_scale = rep(NA_real_, n),
    gp_shape = NA_real_,
    gp_shape_se = NA_real_,
    n_hours_used = 0L,
    n_eff_median = NA_real_
  )
  if (!any(usable)) {
    return(empty)
  }

  ok <- pieces[usable]
  excess <- lapply(ok, `[[`, "excess")
  weight <- lapply(ok, `[[`, "weight")

  if (is.null(shape)) {
    fit <- fit_gpd_shared_shape(
      excess,
      weight,
      n_boot = n_boot,
      bias_correct = bias_correct
    )
    xi <- fit$shape
    xi_se <- fit$shape_se
    scale_ok <- fit$scale
  } else {
    # Shape supplied (e.g. smoothed across neighbours): only the scales are free.
    xi <- shape
    xi_se <- NA_real_
    scale_ok <- vapply(
      seq_along(excess),
      function(i) gpd_profile_scale(excess[[i]], weight[[i]], xi)$scale,
      numeric(1L)
    )
  }
  if (!is.finite(xi)) {
    return(empty)
  }

  out <- empty
  out$graft_of[usable] <- vapply(ok, `[[`, numeric(1L), "threshold")
  out$graft_tail_prob[usable] <- vapply(ok, `[[`, numeric(1L), "tail_prob")
  out$gp_scale[usable] <- scale_ok
  out$gp_shape <- xi
  out$gp_shape_se <- xi_se
  out$n_hours_used <- sum(usable)
  out$n_eff_median <- stats::median(
    vapply(ok, `[[`, numeric(1L), "n_eff"),
    na.rm = TRUE
  )
  out
}

#' Build a mixture tail from encoded peak-hour predictions
#'
#' Takes the compact columns written by [dl_encode_peak_hour_distributions()]
#' and returns the closed-form tail of their equal-weight mixture, without
#' rebuilding a single distribution object.
#'
#' @param df A data frame with `graft_of`, `gp_scale`, `gp_shape` and, where
#'   available, `graft_tail_prob`. When `graft_tail_prob` is absent it defaults
#'   to `adaptive_threshold`, which is what the graft rule targets; supply the
#'   column for exact results.
#' @param adaptive_threshold Fallback tail probability, see above.
#' @param compress Number of representative components to keep, or `NULL` for
#'   none. See [mixture_tail_compress()].
#' @returns A `mixture_tail`.
#' @seealso [mixture_tail_return_level()]
#' @export
dl_cell_mixture_tail <- function(
  df,
  adaptive_threshold = 0.5,
  compress = 64L
) {
  needed <- c("graft_of", "gp_scale", "gp_shape")
  missing <- setdiff(needed, names(df))
  if (length(missing)) {
    stop(
      "Encoded predictions are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  tail_prob <- if ("graft_tail_prob" %in% names(df)) {
    df$graft_tail_prob
  } else {
    rep(1 - adaptive_threshold, nrow(df))
  }

  mt <- mixture_tail(
    threshold = df$graft_of,
    tail_prob = tail_prob,
    scale = df$gp_scale,
    shape = df$gp_shape
  )
  if (!is.null(compress) && isTRUE(mt$shared_shape)) {
    mt <- mixture_tail_compress(mt, compress)
  }
  mt
}
