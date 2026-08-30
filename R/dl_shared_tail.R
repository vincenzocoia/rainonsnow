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

#' Shared tail from standardised observations
#'
#' Estimates one tail shape for a cell by standardising each observed response
#' by its OWN predicted location and scale, pooling the standardised values, and
#' fitting a single generalized Pareto tail to their exceedances.
#'
#' Two things make this different from pooling the predictive distributions
#' themselves ([fit_gpd_shared_shape()]).
#'
#' It removes the incidental-parameter problem. There is no per-hour scale left
#' free in the likelihood: each hour's scale comes from a wide, well-determined
#' gap between two central quantiles of its predictive distribution, which the
#' whole forest informs, rather than from the handful of points in that hour's
#' tail.
#'
#' And it is a real likelihood. A forest's predictive distributions at different
#' hours are re-weightings of the same training responses, so pooling them counts
#' every observation many times over -- the point estimate survives that, but
#' nothing derived from the likelihood's curvature does. The standardised values
#' here are one per observation, so standard errors, likelihood-ratio intervals
#' and likelihood-ratio tests all mean what they usually mean.
#'
#' @section Choosing the quantiles:
#' `p_location` and `p_scale` must straddle a gap wide enough to be well
#' determined but near enough the tail to be relevant to it. In simulation the
#' 70th and 95th percentiles beat both a lower pair (50th/90th) and a higher one
#' (80th/98th); the higher pair is worst, because it makes the normaliser itself
#' a tail estimate, which is the thing being avoided. The POT threshold is set
#' separately on the standardised scale via `threshold_prob`, and raising it
#' costs a lot of variance for no reliable gain in bias.
#'
#' @param distributions List of predictive distributions for one cell.
#' @param observed Numeric vector of the responses actually observed at those
#'   same rows, one per distribution.
#' @param p_location,p_scale Probabilities defining each hour's location and
#'   scale: location is the `p_location` quantile, scale is the gap between the
#'   `p_scale` and `p_location` quantiles.
#' @param threshold_prob Quantile of the pooled standardised values at which the
#'   generalized Pareto tail starts.
#' @returns A list with per-hour `graft_of` and `gp_scale`, and the scalars
#'   `gp_shape`, `gp_shape_se` (from the pooled fit's own likelihood) and
#'   `n_exceedances`.
#' @seealso [dl_fit_cell_shared_tail()]
#' @export
fit_shared_tail_standardised <- function(
  distributions,
  observed,
  p_location = 0.7,
  p_scale = 0.95,
  threshold_prob = 0.6
) {
  n <- length(distributions)
  empty <- list(
    graft_of = rep(NA_real_, n),
    gp_scale = rep(NA_real_, n),
    gp_shape = NA_real_,
    gp_shape_se = NA_real_,
    n_exceedances = 0L
  )
  if (length(observed) != n) {
    stop("`observed` must have one value per distribution.", call. = FALSE)
  }

  q_at <- function(dst, p) {
    if (is.null(dst) || !inherits(dst, "dst")) {
      return(NA_real_)
    }
    v <- try(distionary::eval_quantile(dst, at = p), silent = TRUE)
    if (inherits(v, "try-error") || length(v) != 1L) NA_real_ else as.numeric(v)
  }

  loc <- vapply(distributions, q_at, numeric(1L), p_location)
  hi <- vapply(distributions, q_at, numeric(1L), p_scale)
  sc <- hi - loc

  ok <- is.finite(loc) & is.finite(sc) & sc > 0 & is.finite(observed)
  if (sum(ok) < 20L) {
    return(empty)
  }

  z <- (observed[ok] - loc[ok]) / sc[ok]
  t0 <- stats::quantile(z, threshold_prob, names = FALSE)
  excess <- z[z > t0] - t0
  if (length(excess) < 10L) {
    return(empty)
  }

  fit <- fit_gpd_weighted(excess, rep(1, length(excess)))
  if (!is.finite(fit$shape) || !is.finite(fit$scale)) {
    return(empty)
  }

  # Back to the original scale. y = loc + sc * z, so an excess in z maps to
  # sc * excess in y: the shape is unchanged and the scale is multiplied by sc.
  out <- empty
  out$graft_of[ok] <- loc[ok] + t0 * sc[ok]
  out$gp_scale[ok] <- sc[ok] * fit$scale
  out$gp_shape <- fit$shape
  # A genuine likelihood over independent exceedances, so the usual asymptotic
  # standard error applies: sqrt((1 + xi)^2 / k).
  out$gp_shape_se <- sqrt((1 + fit$shape)^2 / length(excess))
  out$n_exceedances <- length(excess)
  out
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
#' @param shape Optionally fix the shape instead of estimating it, while still
#'   refitting each component's scale. Used to apply one cell's fitted shape to a
#'   grid of covariate values, as `apps/7-rain-snow-given-runoff` does.
#' @param method How to estimate the shared shape. `"standardise"` (recommended)
#'   standardises each observed response by its own predicted location and scale
#'   and fits one generalized Pareto to the pooled exceedances -- see
#'   [fit_shared_tail_standardised()]; it requires `observed`. `"pool_predictive"`
#'   pools the predictive distributions themselves with a free scale each.
#' @param observed Numeric vector of observed responses, one per distribution.
#'   Required by `method = "standardise"`.
#' @param ... Passed to [fit_shared_tail_standardised()].
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
  bias_correct = FALSE,
  method = c("pool_predictive", "standardise"),
  observed = NULL,
  ...
) {
  method <- match.arg(method)

  # The standardised route needs the observed responses, and only estimates the
  # shape; a caller supplying `shape` wants the other branch, which refits each
  # component's scale under a shape it already has.
  if (identical(method, "standardise") && is.null(shape)) {
    if (is.null(observed)) {
      stop(
        "`method = \"standardise\"` needs `observed`, the response actually ",
        "seen at each row.",
        call. = FALSE
      )
    }
    std <- fit_shared_tail_standardised(distributions, observed, ...)
    if (is.finite(std$gp_shape)) {
      keep <- is.finite(std$graft_of) & is.finite(std$gp_scale)
      tp <- rep(NA_real_, length(distributions))
      tp[keep] <- vapply(
        which(keep),
        function(i) {
          v <- try(
            distionary::eval_survival(distributions[[i]], at = std$graft_of[i]),
            silent = TRUE
          )
          if (inherits(v, "try-error") || length(v) != 1L) NA_real_ else as.numeric(v)
        },
        numeric(1L)
      )
      return(list(
        graft_of = std$graft_of,
        graft_tail_prob = tp,
        gp_scale = std$gp_scale,
        gp_shape = std$gp_shape,
        gp_shape_se = std$gp_shape_se,
        n_hours_used = sum(keep),
        n_eff_median = std$n_exceedances
      ))
    }
    # Fall through to the pooled route if the standardised fit failed.
  }

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
    # Shape supplied by the caller: only the scales are free.
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
