# Distributional learning (quantile regression forest per cell) on POT peak hours.
# Model formula: inputs/distributional_learning.yaml
# Writes: derived/era5_land_hourly_alps_dl_{rqforest_models,predictions,diagnostics}.rds
#   (predictions RDS is compact: GPD tail stored as numeric columns, not list-column dsts)
# Downstream: marginal mixtures + return levels — scripts/5-runoff_marginals.r;
# joint rainfall–snowmelt — scripts/6-drivers_joint_distribution.r;
# conditional rain–snow given runoff — scripts/7-likeliest_rain_snow.r
# %%
library(tidyverse)
library(rlang)
library(yaml)
library(probaverse)
library(logger)
devtools::load_all()

meta <- read_yaml(here::here("inputs", "distributional_learning.yaml"))
rq <- meta$dl_rqforest
dl_fit_args <- c(
  list(
    yname = rq$yname,
    xnames = unlist(rq$xnames, use.names = FALSE),
    na_action = rq$na_action %||% "drop",
    min_obs = as.integer(rq$min_obs %||% 5L)
  ),
  rq[setdiff(names(rq), c("yname", "xnames", "na_action", "min_obs"))]
)

log_info("Starting 4-distributional_learning.r")

dat <- read_rds(here::here("derived", "era5_land_hourly_alps_peaks.rds"))

# %%
log_info("Fitting dl_rqforest at each grid cell")
models <- dat |>
  nest(data = !c(cell_id, x, y)) |>
  mutate(
    dl_rqforest = map(
      data,
      \(cell_data) do.call(dl_rqforest, c(list(data = cell_data), dl_fit_args)),
      .progress = TRUE
    )
  )

log_info(
  "Writing per-cell dl_rqforest models (for apps/5-return-level-explorer)"
)
models |>
  select(cell_id, x, y, dl_rqforest) |>
  write_rds(here::here(
    "derived",
    "era5_land_hourly_alps_dl_rqforest_models.rds"
  ))

# %%
log_info(
  "Predictive distributions per hour: rqforest conditional dist and GP-tail version"
)
forest_models <- models |>
  mutate(
    distribution_forest = map2(dl_rqforest, data, predict, .progress = TRUE)
  ) |>
  select(!dl_rqforest) |>
  unnest(c(data, distribution_forest))

# %%
# GP tail per peak hour.
#
# With `shape_pooling: shared` (the default) the shape is fitted once per cell
# and only the scale varies by hour. Fitting a shape per hour instead makes the
# cell's marginal -- an equal-weight mixture of these -- inherit max_i(xi_i),
# which is the maximum of a few hundred noisy estimates. See
# scripts/experiments/tail-index-pooling.R for what that costs.
tail_cfg <- meta$gp_tail %||% list()
shape_pooling <- tail_cfg$shape_pooling %||% "shared"
adaptive_threshold <- as.numeric(tail_cfg$adaptive_threshold %||% 0.5)
n_boot <- as.integer(tail_cfg$n_boot %||% 25L)
bias_correct <- isTRUE(tail_cfg$bias_correct %||% FALSE)
smooth_cfg <- tail_cfg$spatial_smoothing %||% list()
smooth_radius <- as.integer(smooth_cfg$radius %||% 0L)

# The bootstrap only produces the shape's standard error and, optionally, its
# bias correction. The standard error is consumed by the spatial smoothing and
# nothing else, so when smoothing is off and the correction is off there is
# nothing for 25 extra fits per cell to do.
if (!bias_correct && smooth_radius <= 0L) {
  n_boot <- 0L
}

log_info(paste("GP tail: shape_pooling =", shape_pooling))

if (identical(shape_pooling, "shared")) {
  cell_tails <- forest_models |>
    nest(cell_rows = !c(cell_id, x, y)) |>
    mutate(
      tail_fit = map(
        cell_rows,
        \(df) dl_fit_cell_shared_tail(
          df$distribution_forest,
          adaptive_threshold = adaptive_threshold,
          n_boot = n_boot,
          bias_correct = bias_correct
        ),
        .progress = TRUE
      )
    )

  # Borrow shape between neighbouring cells, then refit each hour's scale under
  # the smoothed shape. A cell with a precise estimate barely moves. This is the
  # only step that lets one cell influence another, and it is off by default.
  if (smooth_radius > 0L && nrow(cell_tails) > 2L) {
    raw_shape <- map_dbl(cell_tails$tail_fit, "gp_shape")
    raw_se <- map_dbl(cell_tails$tail_fit, "gp_shape_se")
    smoothed <- smooth_tail_shape(
      raw_shape,
      raw_se,
      cell_tails$x,
      cell_tails$y,
      radius = smooth_radius,
      min_neighbours = as.integer(smooth_cfg$min_neighbours %||% 2L)
    )
    log_info(sprintf(
      "Spatial shape smoothing: between-cell sd = %.4f; mean weight on own estimate = %.2f",
      smoothed$tau,
      mean(smoothed$weight, na.rm = TRUE)
    ))
    cell_tails <- cell_tails |>
      mutate(
        shape_raw = raw_shape,
        shape_smoothed = smoothed$shape,
        tail_fit = map2(
          cell_rows,
          smoothed$shape,
          \(df, xi) dl_fit_cell_shared_tail(
            df$distribution_forest,
            adaptive_threshold = adaptive_threshold,
            shape = xi
          ),
          .progress = TRUE
        )
      )
    write_rds(
      select(cell_tails, cell_id, x, y, shape_raw, shape_smoothed),
      here::here("derived", "era5_land_hourly_alps_dl_tail_shapes.rds")
    )
  }

  peak_hour_distributions <- cell_tails |>
    mutate(
      cell_rows = map2(cell_rows, tail_fit, \(df, tf) {
        mutate(
          df,
          graft_of = tf$graft_of,
          graft_tail_prob = tf$graft_tail_prob,
          gp_scale = tf$gp_scale,
          gp_shape = tf$gp_shape
        )
      })
    ) |>
    select(!tail_fit) |>
    unnest(cell_rows) |>
    mutate(
      # Diagnostics still want the distribution objects; the RDS stays compact.
      distribution_gp = pmap(
        list(distribution_forest, graft_of, gp_scale, gp_shape),
        \(d, u, s, k) {
          if (!is.finite(u) || !is.finite(s) || !is.finite(k)) {
            return(d)
          }
          reconstruct_graft_gp(d, u, s, k)
        },
        .progress = TRUE
      )
    )
} else {
  peak_hour_distributions <- mutate(
    forest_models,
    distribution_gp = map(
      distribution_forest,
      fit_and_graft_gp,
      adaptive_threshold = adaptive_threshold,
      .progress = TRUE
    )
  )
}

peak_hour_distributions <- peak_hour_distributions |>
  select(
    cell_id,
    x,
    y,
    date,
    rainfall_hourly,
    snowmelt_hourly,
    runoff_hourly,
    contains("distribution"),
    any_of(c("graft_of", "graft_tail_prob", "gp_scale", "gp_shape"))
  )

# %%
log_info("Diagnostics (P-P calibration, quantile skill vs marginal)")
dl_diag <- dl_build_diagnostics(peak_hour_distributions, dat)
diag_path <- here::here("derived", dl_diagnostics_basename())
dl_write_diagnostics(dl_diag, diag_path)
log_info(paste("Wrote", diag_path))

plots_dir <- here::here("plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

p_pp <- ggplot(dl_diag$pp, aes(p_empirical, p_model)) +
  facet_wrap(~model, nrow = 1) +
  geom_line(aes(group = cell_id), alpha = 0.1) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    colour = "orange3"
  ) +
  theme_bw()

ggplot2::ggsave(
  file.path(plots_dir, "dl_pp_calibration.pdf"),
  p_pp,
  width = 7.2,
  height = 5.4
)

p_skill <- ggplot(dl_diag$skill, aes(tau, skill_score)) +
  geom_line(
    aes(group = interaction(x, y)),
    alpha = 0.33
  ) +
  labs(
    x = "Quantile Level",
    title = "Skill Score of Quantile Regression Forest Method"
  ) +
  scale_y_continuous("Skill", labels = scales::percent_format()) +
  theme_bw()

ggplot2::ggsave(
  file.path(plots_dir, "dl_quantile_skill_score.pdf"),
  p_skill,
  width = 7.2,
  height = 5.4
)

log_info(paste("Wrote", file.path(plots_dir, "dl_pp_calibration.pdf")))
log_info(paste("Wrote", file.path(plots_dir, "dl_quantile_skill_score.pdf")))

# %%
log_info(
  "Writing hourly DL predictions (compact RDS; GPD tail as numeric columns)"
)
predictions_path <- here::here(
  "derived",
  "era5_land_hourly_alps_dl_predictions.rds"
)
encoded <- dl_encode_peak_hour_distributions(peak_hour_distributions)
log_info(
  paste0(
    "Predictions RDS payload: ",
    format(object.size(encoded), units = "auto"),
    " in memory (",
    nrow(encoded),
    " rows)"
  )
)
dl_write_peak_hour_predictions(encoded, predictions_path)

log_info("Finished 4-distributional_learning.r")
