# Marginal distribution modelling: the equal-weight mixture of a cell's
# peak-hour predictive distributions, and the return-level table used by
# apps/5-runoff-marginals-explorer and script 7.
#
# Requires: outputs of scripts/3-pot_spatial_eo.r and 4-distributional_learning.r
#
# HOW THE MARGINAL IS EVALUATED
#
# Each peak hour's predictive distribution is empirical below its graft point
# and generalized Pareto above it, so the cell mixture splits the same way:
#
#   * Body. A mixture of finite distributions is itself a finite distribution,
#     so mix2() aggregates the weights into a single empirical distribution.
#     Cheap and exact.
#
#   * Tail. Above max(graft_of) every component is in its GP part, and the
#     mixture survival is elementary in the stored (graft_of, graft_tail_prob,
#     gp_scale, gp_shape) columns. dl_cell_mixture_tail() evaluates it as one
#     matrix operation, with no distribution objects built at all -- which is
#     what makes the long return periods fast.
#
# Return periods whose level falls between min(graft_of) and max(graft_of) are
# in a transition band where some components are still empirical; those are
# taken from the body mixture and flagged in the `region` column. Everything at
# the return periods that matter for design (roughly T >= 10 years) sits in the
# exact tail region.
# %%
library(tidyverse)
library(rlang)
library(logger)
library(yaml)
library(probaverse)
devtools::load_all()

# The response the transport step fits its copula against, read from the same
# config the forest was fitted with so the two cannot drift apart.
transport_response <- read_yaml(
  here::here("inputs", "distributional_learning.yaml")
)$dl_rqforest$yname

return_periods <- rp_marginal_curve()

log_info("Starting 5-runoff_marginals.r")

dat <- read_rds(here::here("derived", "era5_land_hourly_alps_peaks.rds"))

# Read WITHOUT decoding: the tail is evaluated from the numeric columns, so
# there is no reason to rebuild grafted distribution objects.
encoded <- dl_read_peak_hour_predictions(
  here::here("derived", "era5_land_hourly_alps_dl_predictions.rds"),
  decode = FALSE
)

if (!dl_predictions_is_encoded(encoded)) {
  log_warn("Predictions file is in the legacy format; encoding in memory.")
  encoded <- dl_encode_peak_hour_distributions(encoded)
}

has_shared_shape <- encoded |>
  group_by(cell_id) |>
  summarise(n_shapes = n_distinct(round(gp_shape, 10)), .groups = "drop") |>
  pull(n_shapes) |>
  max(na.rm = TRUE) |>
  identical(1L)

if (!has_shared_shape) {
  log_warn(paste(
    "Tail shapes vary WITHIN cells (a separate shape per peak hour), so each",
    "cell's marginal inherits max_i(xi_i) -- the noisiest single hour. Set",
    "gp_tail.shape_pooling: shared in inputs/distributional_learning.yaml and",
    "rerun script 4. That pools across a cell's peak hours only; every cell is",
    "still fitted on its own data."
  ))
}

# %%
log_info("Events per year per cell (the POT event axis)")
num_pot_events <- dat |>
  group_by(cell_id, x, y) |>
  summarise(
    num_events_per_year = n() / (diff(range(year(date))) + 1),
    .groups = "drop"
  )

# %%
log_info("Closed-form mixture tails per cell")
cell_tails <- encoded |>
  nest(rows = !c(cell_id, x, y)) |>
  left_join(num_pot_events, by = c("cell_id", "x", "y")) |>
  mutate(
    mixture_tail = map(
      rows,
      \(df) tryCatch(dl_cell_mixture_tail(df), error = function(e) NULL),
      .progress = TRUE
    )
  )

n_failed <- sum(map_lgl(cell_tails$mixture_tail, is.null))
if (n_failed > 0L) {
  log_warn(paste(n_failed, "cells had no usable tail and are dropped."))
}
cell_tails <- filter(cell_tails, !map_lgl(mixture_tail, is.null))

log_info("Writing mixture tails to file")
cell_tails |>
  select(cell_id, x, y, num_events_per_year, mixture_tail) |>
  write_rds(here::here("derived", "era5_land_hourly_alps_dl_mixture_tails.rds"))

# %%
# Body mixture: a mixture of finite distributions collapses to a single
# empirical distribution, so this stays cheap even with many peak hours.
log_info("Body mixtures (empirical) per cell")
if ("distribution_forest" %in% names(encoded)) {
  body_mixtures <- encoded |>
    nest(rows = !c(cell_id, x, y)) |>
    mutate(
      marginal_forest = map(
        rows,
        \(df) tryCatch(
          mix2(df$distribution_forest, na_action_dst = "drop"),
          error = function(e) NULL
        ),
        .progress = TRUE
      )
    ) |>
    select(cell_id, x, y, marginal_forest)

  write_rds(
    body_mixtures,
    here::here("derived", "era5_land_hourly_alps_dl_marginals.rds")
  )
} else {
  log_warn("No distribution_forest column; body mixture unavailable.")
  body_mixtures <- tibble(
    cell_id = integer(),
    x = double(),
    y = double(),
    marginal_forest = list()
  )
}

# %%
log_info("Marginal return levels")

# Two curves per cell, matching the schema apps/5-runoff-marginals-explorer
# expects:
#   "Random Forest"  -- the empirical (body) mixture on its own. It cannot go
#                       past its largest outcome, which is the limitation the GP
#                       tail exists to fix.
#   "GP conversion"  -- the closed-form GP mixture tail, falling back to the
#                       body for return periods below the tail region.
# `region` records which of the two supplied each GP-conversion level.
return_levels_long <- cell_tails |>
  left_join(body_mixtures, by = c("cell_id", "x", "y")) |>
  mutate(
    curve = pmap(
      list(mixture_tail, num_events_per_year, marginal_forest),
      \(mt, nep, body) {
        p <- 1 / (return_periods * nep)
        have_body <- inherits(body, "dst")

        level_body <- rep(NA_real_, length(p))
        if (have_body) {
          # mix2() collapses a mixture of finite distributions to a single
          # empirical one, whose quantile works. If it did not collapse, the
          # result is a generic "Mixture" whose quantile is disabled in the
          # installed distionary -- invert its CDF numerically instead.
          level_body <- tryCatch(
            distionary::eval_quantile(body, at = 1 - p),
            error = function(e) eval_return_numeric(body, at = 1 / p)
          )
        }

        level_tail <- mixture_tail_quantile(mt, p)
        need_body <- is.na(level_tail)

        dplyr::bind_rows(
          tibble(
            return_period = return_periods,
            model = "Random Forest",
            return_level = level_body,
            region = "body"
          ),
          tibble(
            return_period = return_periods,
            model = "GP conversion",
            return_level = coalesce(level_tail, level_body),
            region = if_else(need_body, "body", "tail")
          )
        )
      },
      .progress = TRUE
    )
  ) |>
  select(cell_id, x, y, num_events_per_year, curve) |>
  unnest(curve)

gp_rows <- filter(return_levels_long, model == "GP conversion")
log_info(sprintf(
  "GP-conversion levels: %.1f%% from the closed-form tail, %.1f%% from the body mixture",
  100 * mean(gp_rows$region == "tail"),
  100 * mean(gp_rows$region == "body")
))

log_info("Writing marginal return levels to file")
write_rds(
  return_levels_long,
  here::here("derived", "era5_land_hourly_alps_dl_return_levels.rds")
)

# %%
# Per-cell tail summary: the single GP that matches each cell's mixture in the
# limit. Useful for mapping the tail index and spotting cells whose shape is
# out of line with their neighbours.
log_info("Writing per-cell tail summaries")
cell_tails |>
  mutate(summary = map(mixture_tail, mixture_tail_gpd_approx)) |>
  transmute(
    cell_id,
    x,
    y,
    num_events_per_year,
    tail_threshold = map_dbl(summary, "threshold"),
    tail_prob = map_dbl(summary, "tail_prob"),
    tail_scale = map_dbl(summary, "scale"),
    tail_shape = map_dbl(summary, "shape")
  ) |>
  write_rds(
    here::here("derived", "era5_land_hourly_alps_dl_tail_summary.rds")
  )

# %%
# Copula transport of the marginal, as an alternative to the mixture.
#
# Off unless inputs/distributional_learning.yaml turns it on, and written to its
# own files either way, so nothing above changes shape. The mixture is the law
# of total probability over the covariate values in the sample and therefore
# carries the conditional tail index; the transport pushes each hour's tail back
# through the copula so that every hour estimates the whole marginal, and takes
# the pointwise median of those estimates. See R/dl_transport_marginal.R for the
# one approximation it rests on.
transport_cfg <- read_yaml(
  here::here("inputs", "distributional_learning.yaml")
)$transport

if (isTRUE(transport_cfg$enabled)) {
  family <- transport_cfg$family %||% "survival_clayton"
  combine <- transport_cfg$combine %||% "median"
  spread_at <- as.numeric(transport_cfg$spread_at %||% 1e-4)
  log_info(paste("Copula transport of the marginal: family =", family))

  transported <- encoded |>
    nest(rows = !c(cell_id, x, y)) |>
    left_join(num_pot_events, by = c("cell_id", "x", "y")) |>
    mutate(
      ensemble = map(
        rows,
        \(df) tryCatch(
          dl_cell_transport_ensemble(
            df,
            observed = df[[transport_response]],
            family = family
          ),
          error = function(e) NULL
        ),
        .progress = TRUE
      )
    ) |>
    filter(!map_lgl(ensemble, is.null))

  log_info(sprintf(
    "%d of %d cells produced a transported marginal",
    nrow(transported),
    n_distinct(encoded$cell_id)
  ))

  transported |>
    mutate(
      copula_par = map_dbl(ensemble, "par"),
      n_hours = map_int(ensemble, \(me) length(me$u)),
      transport_spread = map_dbl(ensemble, dl_transport_spread, p = spread_at)
    ) |>
    select(cell_id, x, y, num_events_per_year, copula_par, n_hours,
           transport_spread, ensemble) |>
    write_rds(here::here(
      "derived",
      "era5_land_hourly_alps_dl_transport_marginals.rds"
    ))

  transported |>
    mutate(
      curve = map2(
        ensemble,
        num_events_per_year,
        \(me, nep) dl_transport_return_levels(
          me,
          return_periods,
          events_per_year = nep,
          combine = combine
        )
      )
    ) |>
    select(cell_id, x, y, num_events_per_year, curve) |>
    unnest(curve) |>
    mutate(model = "Copula transport", region = "tail") |>
    write_rds(here::here(
      "derived",
      "era5_land_hourly_alps_dl_transport_return_levels.rds"
    ))

  log_info("Wrote transported marginals and return levels")
} else {
  log_info("Copula transport is off (inputs/distributional_learning.yaml)")
}

# %%
log_info("Finished 5-runoff_marginals.r")
