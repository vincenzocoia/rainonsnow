# %%
library(tidyverse)
library(probaverse)
library(terra)
library(tidyterra)
devtools::load_all()

# "Vibrant Summer" palette https://coolors.co/palette/ff595e-ffca3a-8ac926-1982c4-6a4c93
pal <- rev(c("#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93"))

# %%
# Load netcdf file
ds <- rast(here::here("data", "era5_land_hourly_alps.nc"))

res(ds)
ncell(ds)

# Turn into a tibble
ds_tbl <- as_tibble(ds, xy = TRUE)

# The NetCDF stores latitude in `x` and longitude in `y`.
map_xlim <- range(ds_tbl$y, na.rm = TRUE) + c(-0.5, 0.5)
map_ylim <- range(ds_tbl$x, na.rm = TRUE) + c(-0.5, 0.5)

map_bbox <- sf::st_bbox(
  c(
    xmin = map_xlim[1],
    xmax = map_xlim[2],
    ymin = map_ylim[1],
    ymax = map_ylim[2]
  ),
  crs = sf::st_crs(4326)
)

world_map <- rnaturalearth::ne_countries(
  scale = 50,
  returnclass = "sf"
) |>
  sf::st_crop(map_bbox)

spatial_basemap <- list(
  geom_sf(
    data = world_map,
    inherit.aes = FALSE,
    fill = "grey96",
    color = "grey70",
    linewidth = 0.2
  ),
  coord_sf(
    xlim = map_xlim,
    ylim = map_ylim,
    expand = FALSE
  ),
  labs(x = "Longitude", y = "Latitude"),
  theme_minimal()
)

spatial_tile_basemap <- list(
  ggspatial::annotation_map_tile(
    type = "osm",
    zoomin = -1,
    progress = "none"
  ),
  coord_sf(
    xlim = map_xlim,
    ylim = map_ylim,
    expand = FALSE
  ),
  labs(x = "Longitude", y = "Latitude"),
  theme_minimal()
)

# %%
# Pivot the tibble so that time is its own column
dat <- ds_tbl |> 
  pivot_longer(
    !c(x, y),
    names_to = c(".value", "hour"),
    names_pattern = "^(.+)_([0-9]+)$"
  ) |> 
  mutate(hour = as.numeric(hour))

# Remove near-zero and negatives
dat <- mutate(
  dat,
  total_precipitation_hourly = if_else(total_precipitation_hourly < 1e-7, 0, total_precipitation_hourly),
  snowmelt_hourly = if_else(snowmelt_hourly < 1e-7, 0, snowmelt_hourly)
)

# %%
# Some plots
dat |>
  filter(hour == 14) |> 
  ggplot(aes(y, x)) +
  spatial_basemap +
  geom_tile(aes(fill = runoff_hourly), alpha = 0.85) +
  scale_fill_gradientn(colours = pal)

dat |>
  filter(hour == 14) |> 
  ggplot(aes(y, x)) +
  spatial_tile_basemap +
  geom_tile(aes(fill = runoff_hourly), alpha = 0.65) +
  scale_fill_gradientn(colours = pal)

dat |> 
  filter(x < 46, y > 15) |> 
  ggplot(aes(hour, runoff_hourly)) +
  facet_grid(y ~ x) +
  geom_line()

dat |>
  filter(x < 46, y > 15) |> 
  ggplot(aes(snowmelt_hourly, runoff_hourly)) +
  facet_grid(y ~ x) +
  geom_point(alpha = 0.1)

dat |>
  filter(x < 46, y > 15) |> 
  ggplot(aes(total_precipitation_hourly, runoff_hourly)) +
  facet_grid(y ~ x) +
  geom_point(alpha = 0.1)

# %%
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

# %%
# Predictions on the training data
pred <- models |> 
  mutate(hat = map2(dl_rqforest, data, predict, .progress = TRUE)) |> 
  select(!dl_rqforest) |> 
  unnest(c(data, hat))

# %%
# Take a look at one of the predictive distributions (cdf)
plot(pred$hat[[1]], n = 1000)

# Make a P-P plot for each cell to evaluate model fit.
pp_data <- pred |>
  group_by(x, y) |> 
  mutate(
    p_model = map2(hat, data, \(dst, df) eval_cdf(dst, at = df[["runoff_hourly"]])),
    p_empirical = uscore(p_model)
  ) |> 
  select(x, y, starts_with("p_"))

# %%
# Marginal (Unconditional) distributions
marginals <- pred |> 
  group_by(x, y) |> 
  summarise(marginal_runoff = list(mix(hat)), .groups = "drop")

# Example of one unconditional distribution
plot(marginals$marginal_runoff[[4]], to = 0.00001, n = 1000)
# %%
# Map of 1% and 99% percentile of those unconditional distributions
pctl <- marginals |> 
  mutate(
    df = map(
      marginal_runoff,
      enframe_quantile,
      at = c(0.01, 0.99),
      arg_name = "tau"
    )
  ) |> 
  select(!marginal_runoff) |> 
  unnest(df)

ggplot(pctl, aes(y, x)) +
  spatial_basemap +
  facet_wrap(~ tau, ncol = 2) +
  geom_tile(aes(fill = quantile), alpha = 0.85) +
  scale_fill_gradientn("Percentile:\nHourly Runoff (mm)", colours = pal)

# Map of percentiles on hour 14 vs actual
pctl_hr14 <- pred |> 
  filter(hour == 14) |> 
  mutate(
    runoff_01 = map_dbl(hat, eval_quantile, at = 0.01),
    runoff_99 = map_dbl(hat, eval_quantile, at = 0.99),
    runoff_mean = map_dbl(hat, mean),
  ) |> 
  select(x, y, starts_with("runoff")) |> 
  pivot_longer(
    starts_with("runoff"),
    names_to = c(".value", "quantity"),
    names_sep = "_"
  ) |>
  mutate(quantity = case_when(
    quantity == "01" ~ "1% Percentile",
    quantity == "99" ~ "99% Percentile",
    quantity == "mean" ~ "Mean",
    TRUE ~ "Observed"
  ))

ggplot(pctl_hr14, aes(y, x)) +
  spatial_basemap +
  facet_wrap(~ quantity, nrow = 2) +
  geom_tile(aes(fill = runoff), alpha = 0.85) +
  scale_fill_gradientn("Hourly Runoff (mm)", colours = pal)

# %%
# Percent due to rain
foo <- dat |>
  mutate(
    total_precipitation_hourly = if_else(abs(total_precipitation_hourly) < 1e-7, 0, total_precipitation_hourly),
    snowmelt_hourly = if_else(abs(snowmelt_hourly) < 1e-7, 0, snowmelt_hourly),
    sm = total_precipitation_hourly + snowmelt_hourly,
    pct_precip = total_precipitation_hourly / sm,
    .keep = "used"
  ) |> 
  drop_na()

hist(foo$pct_precip)
