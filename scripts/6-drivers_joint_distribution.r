# Joint distributional modelling of hourly rainfall and snowmelt per cell:
# famish marginals + rvinecopulib bicop. Fitting options: inputs/rain_snow_joint_model.yaml
# Requires: derived/era5_land_hourly_alps_all.rds from scripts/2-tablify_spatial_eo.r
# Downstream: conditional rain–snow given runoff — scripts/7-likeliest_rain_snow.r
# %%
library(tidyverse)
library(rlang)
library(yaml)
library(rvinecopulib)
library(logger)
devtools::load_all()

meta <- read_yaml(here::here("inputs", "rain_snow_joint_model.yaml"))
cfg <- meta$fit_joint_rain_snow_cells

log_info("Starting 6-drivers_joint_distribution.r")

peaks <- read_rds(here::here("derived", "era5_land_hourly_alps_peaks.rds"))

log_info(
  paste(
    "Joint hourly rainfall–snowmelt per cell:",
    "famish::fit_dst() marginals + rvinecopulib::bicop()"
  )
)

# %%
# Single-cell exploration. The focus cell comes from
# `likeliest_rain_snow$cell_id` in inputs/rain_snow_joint_model.yaml; when that is
# null/absent or not present in the current grid, auto-select the cell with the
# most co-occurring rain+snowmelt peak hours (the "mixed regime" cell).
focus_cell <- suppressWarnings(as.integer(meta$likeliest_rain_snow$cell_id %||% NA_integer_))
available_cells <- sort(unique(peaks$cell_id))
if (is.na(focus_cell) || !(focus_cell %in% available_cells)) {
  mixed_counts <- peaks |>
    filter(rainfall_hourly > 0, snowmelt_hourly > 0) |>
    count(cell_id, sort = TRUE)
  focus_cell <- if (nrow(mixed_counts) > 0) {
    mixed_counts$cell_id[[1]]
  } else {
    available_cells[[1]]
  }
  log_info(paste("Focus cell auto-selected (most mixed rain+snow peaks):", focus_cell))
} else {
  log_info(paste("Focus cell from config:", focus_cell))
}

peaks_focus <- filter(peaks, cell_id == focus_cell)

ggplot(peaks_focus, aes(rainfall_hourly, snowmelt_hourly)) +
  geom_point() +
  geom_point(
    data = filter(peaks_focus, rainfall_hourly == 0 | snowmelt_hourly == 0),
    colour = "green",
  ) +
  theme_bw()

fit_joint_rain_snow_cell(
  peaks_focus$rainfall_hourly,
  peaks_focus$snowmelt_hourly
)

# %%
joint_by_cell <- fit_joint_rain_snow_cells(
  peaks,
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


# %%
write_rds(
  joint_by_cell,
  here::here("derived", "era5_land_hourly_alps_joint_rain_snow.rds")
)
log_info("Wrote era5_land_hourly_alps_joint_rain_snow.rds")

log_info("Finished 6-drivers_joint_distribution.r")
