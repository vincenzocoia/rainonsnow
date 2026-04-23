# Distributional learning (quantile regression forest per cell) on POT peak hours.
# Model formula: inputs/dl_metadata.yaml
# Downstream: marginal mixtures + return levels — scripts/5-runoff_marginals.r;
# joint rainfall–snowmelt — scripts/6-drivers_joint_distribution.r;
# conditional rain–snow given runoff — scripts/7-likeliest_rain_snow.r
# %%
library(tidyverse)
library(yaml)
library(probaverse)
library(logger)
devtools::load_all()

meta <- read_yaml(here::here("inputs", "dl_metadata.yaml"))
rq <- meta$dl_rqforest

log_info("Starting 4-distributional_learning.r")

dat <- read_rds(here::here("data", "era5_land_hourly_alps_peaks.rds"))

# %%
log_info("Fitting dl_rqforest at each grid cell")
models <- dat |>
  nest(data = !c(cell_id, x, y)) |>
  mutate(
    dl_rqforest = map(
      data,
      dl_rqforest,
      yname = rq$yname,
      xnames = rq$xnames,
      .progress = TRUE
    )
  )

log_info("Writing per-cell dl_rqforest models (for apps/return-level-explorer)")
models |>
  select(cell_id, x, y, dl_rqforest) |>
  write_rds(here::here("data", "era5_land_hourly_alps_dl_rqforest_models.rds"))

# %%
log_info(
  "Predictive distributions per hour: rqforest conditional dist and GP-tail version"
)
peak_hour_distributions <- models |>
  mutate(
    distribution_forest = map2(dl_rqforest, data, predict, .progress = TRUE)
  ) |>
  select(!dl_rqforest) |>
  unnest(c(data, distribution_forest)) |>
  mutate(
    distribution_gp = map(
      distribution_forest,
      convert_emp_to_gp,
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
  here::here("data", "era5_land_hourly_alps_dl_predictions.rds")
)

# %%
log_info("Diagnostic plots (P-P calibration, quantile skill vs marginal)")
plots_dir <- here::here("plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

pp_data <- peak_hour_distributions |>
  group_by(cell_id, x, y) |>
  mutate(
    p_model_forest = map2_dbl(distribution_forest, runoff_hourly, eval_cdf),
    p_empirical_forest = uscore(p_model_forest),
    p_model_gp = map2_dbl(distribution_gp, runoff_hourly, eval_cdf),
    p_empirical_gp = uscore(p_model_gp),
  ) |>
  ungroup() |>
  select(cell_id, x, y, starts_with("p_")) |>
  pivot_longer(
    starts_with("p_"),
    names_to = c(".value", "model"),
    names_pattern = "(p_.*)_(.*)"
  )

p_pp <- ggplot(pp_data, aes(p_empirical, p_model)) +
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

# %%
qscores_model <- peak_hour_distributions |>
  mutate(
    df = map(
      distribution_forest,
      enframe_quantile,
      at = 1:99 / 100,
      arg_name = "tau"
    )
  ) |>
  unnest(df) |>
  group_by(cell_id, x, y, tau) |>
  summarise(
    qscore_model = mean(quantile_score(
      runoff_hourly,
      xhat = quantile,
      tau = tau
    )),
    .groups = "drop"
  )

null_model <- dat |>
  group_by(cell_id, x, y) |>
  summarise(runoff_hourly = list(runoff_hourly), .groups = "drop") |>
  mutate(marginal = map(runoff_hourly, dst_empirical))

null_quantiles <- null_model |>
  mutate(
    df = map(marginal, enframe_quantile, at = 1:99 / 100, arg_name = "tau")
  ) |>
  select(!marginal) |>
  unnest(df)

qscores_null <- null_quantiles |>
  mutate(
    qscore_null = pmap_dbl(
      list(runoff_hourly, quantile, tau),
      \(y, q, p) mean(quantile_score(y, xhat = q, tau = p))
    )
  ) |>
  select(x, y, tau, qscore_null)

qscores <- left_join(qscores_null, qscores_model, by = c("x", "y", "tau")) |>
  mutate(skill_score = 1 - qscore_model / qscore_null)

p_skill <- ggplot(qscores, aes(tau, skill_score)) +
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
