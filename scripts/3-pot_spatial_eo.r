# Peaks-over-threshold extraction from tabulated EO hourly data
# (output of 2-tablify_spatial_eo.r).
# %%
library(tidyverse)
library(logger)
library(yaml)
devtools::load_all()
log_info("Starting 3-pot_spatial_eo.r")

meta_path <- here::here("inputs", "pot_metadata.yaml")
if (!file.exists(meta_path)) {
  meta <- list(quantile = 0.995, min_gap = 6)
  write_yaml(meta, meta_path)
  log_info(paste("Created default POT metadata at", meta_path))
} else {
  meta <- read_yaml(meta_path)
}

quantile_level <- if (is.null(meta$quantile)) 0.995 else as.numeric(meta$quantile)
min_gap_obs <- if (is.null(meta$min_gap)) {
  6L
} else {
  as.integer(max(0, round(as.numeric(meta$min_gap))))
}
if (length(quantile_level) != 1L || !is.finite(quantile_level)) {
  stop("`quantile` in inputs/pot_metadata.yaml must be a single finite number", call. = FALSE)
}
if (quantile_level <= 0 || quantile_level >= 1) {
  stop("`quantile` in inputs/pot_metadata.yaml must be strictly between 0 and 1", call. = FALSE)
}
if (length(min_gap_obs) != 1L || is.na(min_gap_obs) || min_gap_obs < 0) {
  stop("`min_gap` in inputs/pot_metadata.yaml must be a non-negative integer", call. = FALSE)
}

log_info(paste(
  "POT settings: threshold quantile =", quantile_level,
  "; min_gap =", min_gap_obs, "observations"
))

dat <- read_rds(here::here("data", "era5_land_hourly_alps_all.rds"))

threshold <- function(x) unname(quantile(x, quantile_level, na.rm = TRUE))

# %%
log_info("Taking peaks over thresholds")
pot_nested <- dat |>
  nest(data = !c(cell_id, x, y)) |>
  mutate(
    data = map(
      data,
      \(df) {
        get_pot_events(
          df,
          threshold = threshold,
          flow_col = "runoff_hourly",
          date_col = "date",
          min_gap = min_gap_obs
        )
      },
      .progress = TRUE
    ),
    threshold = map_dbl(data, \(tbl) attr(tbl, "threshold"))
  )

pot_thresholds <- select(pot_nested, cell_id, x, y, threshold)

pot <- unnest(pot_nested, data)

# %%
log_info("Writing peaks to CSV and RDS")
write_csv(pot, file = here::here("data", "era5_land_hourly_alps_peaks.csv"))
write_rds(pot, file = here::here("data", "era5_land_hourly_alps_peaks.rds"))

log_info("Writing POT thresholds (one row per cell)")
write_csv(
  pot_thresholds,
  file = here::here("data", "era5_land_hourly_alps_pot_thresholds.csv")
)
write_rds(
  pot_thresholds,
  file = here::here("data", "era5_land_hourly_alps_pot_thresholds.rds")
)
