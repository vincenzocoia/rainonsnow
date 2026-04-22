# Joint distributional modelling of hourly rainfall and snowmelt per cell:
# famish marginals + rvinecopulib bicop. Fitting options: inputs/joint_rain_snow_metadata.yaml
# Requires: data/era5_land_hourly_alps_all.rds from scripts/2-tablify_spatial_eo.r
# Downstream: conditional rain–snow given runoff — scripts/7-rain_snow_conditional_runoff_spatial_eo.r
# %%
library(tidyverse)
library(rlang)
library(yaml)
library(rvinecopulib)
library(logger)
devtools::load_all()

meta <- read_yaml(here::here("inputs", "joint_rain_snow_metadata.yaml"))
cfg <- meta$fit_joint_rain_snow_cells

log_info("Starting 6-joint_rain_snow_spatial_eo.r")

hourly_all <- read_rds(here::here("data", "era5_land_hourly_alps_all.rds"))

log_info(
  paste(
    "Joint hourly rainfall–snowmelt per cell:",
    "famish::fit_dst() marginals + rvinecopulib::bicop()"
  )
)

joint_by_cell <- fit_joint_rain_snow_cells(
  hourly_all,
  group_cols = cfg$group_cols,
  rainfall_col = cfg$rainfall_col,
  snowmelt_col = cfg$snowmelt_col,
  marginal_rainfall = cfg$marginal_rainfall,
  marginal_snowmelt = cfg$marginal_snowmelt,
  bicop_family_set = cfg$bicop_family_set,
  bicop_controls = cfg$bicop_controls %||% list(),
  min_obs = as.integer(cfg$min_obs %||% 40L),
  progress = isTRUE(cfg$progress),
  verbose = isTRUE(cfg$verbose %||% FALSE)
)

write_rds(
  joint_by_cell,
  here::here("data", "era5_land_hourly_alps_joint_rain_snow.rds")
)
log_info("Wrote era5_land_hourly_alps_joint_rain_snow.rds")

log_info("Finished 6-joint_rain_snow_spatial_eo.r")
