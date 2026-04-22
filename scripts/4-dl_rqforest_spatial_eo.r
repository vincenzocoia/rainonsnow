# Distributional learning (quantile regression forest per cell) on POT peak hours.
# Model formula: inputs/dl_metadata.yaml
# Downstream: marginal mixtures + return levels — scripts/5-marginal_mixtures_spatial_eo.r;
# joint rainfall–snowmelt — scripts/6-joint_rain_snow_spatial_eo.r;
# conditional rain–snow given runoff — scripts/7-rain_snow_conditional_runoff_spatial_eo.r
# %%
library(tidyverse)
library(yaml)
library(probaverse)
library(logger)
devtools::load_all()

meta <- read_yaml(here::here("inputs", "dl_metadata.yaml"))
rq <- meta$dl_rqforest

log_info("Starting 4-dl_rqforest_spatial_eo.r")

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

log_info("Finished 4-dl_rqforest_spatial_eo.r")
