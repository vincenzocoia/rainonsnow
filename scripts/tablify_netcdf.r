library(tidyverse)
library(terra)
library(tidyterra)

ds <- rast(here::here("data", "era5_land_hourly_alps.nc"))

dat <- as_tibble(ds, xy = TRUE) |>
  pivot_longer(
    !c(x, y),
    names_to = c(".value", "hour"),
    names_pattern = "^(.+)_([0-9]+)$"
  ) |>
  mutate(hour = as.numeric(hour))

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

write_csv(dat, file = here::here("data", "era5_land_hourly_alps.csv"))
write_rds(dat, file = here::here("data", "era5_land_hourly_alps.rds"))
