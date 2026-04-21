# Tablify netCDF spatial Earth Observation data for the Alps bbox.
# %%
library(tidyverse)
library(terra)
library(tidyterra)
library(tidyverse)
library(logger)
library(fs)
log_info("Starting 2-tablify_spatial_eo.r")
epsilon <- 1e-7

files <- dir_ls(here::here("data", "eo"), glob = "*.nc")

tmp_dir <- here::here("data", "intermediate")
dir_create(tmp_dir)

# %%
# Read in each file and turn into a tibble
log_info("Reading in each file and turning into a tibble")

for (file in files) {
  year <- as.numeric(unname(str_extract(file, pattern = "[0-9]{4}")))
  log_info(paste("Reading", file))
  t0 <- Sys.time()
  ds <- rast(file)
  this_dat <- as_tibble(ds, xy = TRUE) |>
    pivot_longer(
      !c(x, y),
      names_to = c(".value", "hour"),
      names_pattern = "^(.+)_([0-9]+)$"
    ) |>
    mutate(
      hour = as.numeric(hour),
      date = ymd_hms(paste0(.env$year, "-01-01 00:00:00")) + dhours(hour),
      .after = "y"
    ) |>
    select(!hour)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  write_rds(this_dat, file = path(tmp_dir, paste0(year, ".rds")))
  log_info(paste0("Finished ", basename(file), " in ", round(secs, 2), " s"))
}

# %%
# Bind the tibbles together
log_info("Binding the tibbles together")
tmp_files <- dir_ls(tmp_dir)

dats <- map(tmp_files, read_rds, .progress = TRUE)

dat <- bind_rows(dats)

# %%
# Label cells (check with plot)
xy_combos <- dat |>
  distinct(x, y) |>
  arrange(desc(x), y) |>
  mutate(cell_id = 1:n())

ggplot(xy_combos, aes(y, x)) +
  geom_label(aes(label = cell_id)) +
  theme_bw()

dat <- dat |>
  left_join(xy_combos, by = c("x", "y")) |>
  select(cell_id, everything())

# %% Calculate Rainfall
dat <- mutate(
  dat,
  rainfall_hourly = total_precipitation_hourly - snowfall_hourly
)

# %%
# Remove near-zero and negatives. Units to mm.
log_info("Removing near-zero and negatives")
dat <- mutate(
  dat,
  rainfall_hourly = if_else(
    rainfall_hourly < epsilon,
    0,
    rainfall_hourly * 1000
  ),
  snowmelt_hourly = if_else(
    snowmelt_hourly < epsilon,
    0,
    snowmelt_hourly * 1000
  )
)

# %%
# Write all data to CSV and RDS
log_info("Writing all data to CSV and RDS")
write_csv(dat, file = here::here("data", "era5_land_hourly_alps_all.csv"))
write_rds(dat, file = here::here("data", "era5_land_hourly_alps_all.rds"))
