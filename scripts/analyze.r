# %%
library(tidyverse)
library(probaverse)
library(terra)
library(tidyterra)
devtools::load_all()

# %%
# Load netcdf file
ds <- rast(here::here("data", "era5_land_hourly_alps.nc"))

res(ds)
ncell(ds)

# Turn into a tibble
ds_tbl <- as_tibble(ds, xy = TRUE)

# %%
# Pivot the tibble so that time is its own column
dat <- ds_tbl |> 
  pivot_longer(
    !c(x, y),
    names_to = c(".value", "hour"),
    names_pattern = "^(.+)_([0-9]+)$"
  ) |> 
  mutate(hour = as.numeric(hour))

# Some plots
dat |>
  filter(hour == 14) |> 
  ggplot(aes(x, y)) +
  geom_tile(aes(fill = runoff_hourly))

dat |> 
  filter(x < 46, y > 15) |> 
  ggplot(aes(hour, runoff_hourly)) +
  facet_grid(y ~ x) +
  geom_line()

dat |>
  #group_by(cell) |> 
  #mutate(across(!hour, nscore)) |> 
  filter(x < 46, y > 15) |> 
  ggplot(aes(snowmelt_hourly, runoff_hourly)) +
  facet_grid(y ~ x) +
  geom_point(alpha = 0.1)

dat |>
  #group_by(cell) |> 
  #mutate(across(!hour, nscore)) |> 
  filter(x < 46, y > 15) |> 
  ggplot(aes(total_precipitation_hourly, runoff_hourly)) +
  facet_grid(y ~ x) +
  geom_point(alpha = 0.1)

# Random Forest distributional learning method
models <- dat |>
  nest(data = !c(x, y)) |>
  mutate(dl_rqforest = map(
    data,
    dl_rqforest,
    yname = "runoff_hourly",
    xnames = c("total_precipitation_hourly", "snowmelt_hourly"),
    .progress = TRUE
  ))

pred <- mutate(models, hat = map2(dl_rqforest, data, predict, .progress = TRUE))
