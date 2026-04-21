# %%
library(tidyverse)
library(terra)
library(tidyterra)
library(tidyverse)
library(logger)
library(fs)
devtools::load_all()
log_info("Starting tablify_netcdf.r")

files <- dir_ls(here::here("data"), glob = "*.nc")

# %%
# Read in each file and turn into a tibble
log_info("Reading in each file and turning into a tibble")
dats <- map(
  files[1:2],
  function(file) {
    ds <- rast(file)
    as_tibble(ds, xy = TRUE) |>
      pivot_longer(
        !c(x, y),
        names_to = c(".value", "hour"),
        names_pattern = "^(.+)_([0-9]+)$"
      ) |>
      mutate(hour = as.numeric(hour))
  },
  .progress = TRUE
)

# %%
# Bind the tibbles together
log_info("Binding the tibbles together")
dat <- bind_rows(dats)

# %%
# Remove near-zero and negatives
log_info("Removing near-zero and negatives")
dat <- mutate(
  dat,
  total_precipitation_hourly = if_else(
    total_precipitation_hourly < 1e-7,
    0,
    total_precipitation_hourly
  ),
  snowmelt_hourly = if_else(
    snowmelt_hourly < 1e-7,
    0,
    snowmelt_hourly
  )
)

# %%
# Write all data to CSV and RDS
log_info("Writing all data to CSV and RDS")
write_csv(dat, file = here::here("data", "era5_land_hourly_alps_all.csv"))
write_rds(dat, file = here::here("data", "era5_land_hourly_alps_all.rds"))


# %%
# Take peaks over thresholds
log_info("Taking peaks over thresholds")
pot <- get_pot_events(
  dat,
  threshold = function(x) unname(quantile(x, 0.8, na.rm = TRUE)),
  flow_col = "runoff_hourly",
  date_col = "date"
)

# Write peaks to CSV and RDS
log_info("Writing peaks to CSV and RDS")
write_csv(pot, file = here::here("data", "era5_land_hourly_alps_peaks.csv"))
write_rds(pot, file = here::here("data", "era5_land_hourly_alps_peaks.rds"))
