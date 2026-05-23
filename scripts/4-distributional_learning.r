# Distributional learning (quantile regression forest per cell) on POT peak hours.
# Model formula: inputs/distributional_learning.yaml
# Writes: derived/era5_land_hourly_alps_dl_{rqforest_models,predictions,diagnostics}.rds
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

peak_hour_distributions <- mutate(
  forest_models,
  distribution_gp = map(
    distribution_forest,
    fit_and_graft_gp,
    adaptive_threshold = 0.5,
    .progress = TRUE
  )
)

peak_hour_distributions <- peak_hour_distributions |>
  select(
    cell_id,
    x,
    y,
    date,
    rainfall_hourly,
    snowmelt_hourly,
    runoff_hourly,
    contains("distribution")
  )

# %%
log_info("Writing hourly DL predictions (single RDS)")
write_rds(
  peak_hour_distributions,
  here::here("derived", "era5_land_hourly_alps_dl_predictions.rds")
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

log_info("Finished 4-distributional_learning.r")
