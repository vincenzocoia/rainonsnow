# Peaks-over-threshold extraction from tabulated EO hourly data
# (output of 2-tablify_spatial_eo.r).
# %%
library(tidyverse)
library(logger)
devtools::load_all()
log_info("Starting 3-pot_spatial_eo.r")

dat <- read_rds(here::here("data", "era5_land_hourly_alps_all.rds"))

threshold <- function(x) unname(quantile(x, 0.995, na.rm = TRUE))

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
          min_gap = 6
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
