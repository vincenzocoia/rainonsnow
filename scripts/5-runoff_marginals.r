# Marginal distribution modelling: equal-weight mixture of peak-hour predictive
# distributions per cell (forest + GP tail) and precomputed return-level table
# for apps/return-level-explorer.
# Requires: outputs of scripts/3-pot_spatial_eo.r and 4-distributional_learning.r
# %%
library(tidyverse)
library(logger)
library(probaverse)
devtools::load_all()

return_periods <- rp_reporting()

log_info("Starting 5-runoff_marginals.r")

dat <- read_rds(here::here("data", "era5_land_hourly_alps_peaks.rds"))
peak_hour_distributions <- read_rds(
  here::here("data", "era5_land_hourly_alps_dl_predictions.rds")
)

# %%
log_info("Marginal mixtures (forest + GP) per cell via mix2()")
marginals_tbl <- peak_hour_distributions |>
  nest(data = !c(cell_id, x, y)) |>
  mutate(
    marginal_forest = map(
      data,
      \(df) {
        mix2(df$distribution_forest, na_action_dst = "drop")
      },
      .progress = TRUE
    ),
    marginal_gp = map(
      data,
      \(df) {
        mix2(df$distribution_gp, na_action_dst = "drop")
      },
      .progress = TRUE
    )
  )

write_rds(
  marginals_tbl,
  here::here("data", "era5_land_hourly_alps_dl_marginals.rds")
)

# %%
log_info("Precomputing marginal return levels (matches Shiny slow-path logic)")
num_pot_events <- dat |>
  group_by(cell_id, x, y) |>
  summarise(
    num_events_per_year = n() / (diff(range(year(date))) + 1),
    .groups = "drop"
  )

marginals_tbl <- marginals_tbl |>
  left_join(num_pot_events, by = c("cell_id", "x", "y"))

# Marginals are evaluated on an event-frequency axis; see enframe_at_events().

return_levels <- mutate(
  marginals_tbl,
  return_period = list(.env$return_periods),
  levels_forest = map2(
    marginal_forest,
    num_events_per_year,
    \(dist, num) eval_return(dist, at = return_periods * num),
    .progress = TRUE
  )
)

return_levels <- mutate(
  return_levels,
  levels_gp = map2(
    marginal_gp,
    num_events_per_year,
    \(dist, num) eval_return(dist, at = return_periods * num),
    .progress = TRUE
  )
)

return_levels <- return_levels |>
  select(!c(starts_with("marginal_"), data)) |>
  unnest(c(return_period, starts_with("levels_")))

return_levels_long <- return_levels |>
  pivot_longer(
    c(levels_forest, levels_gp),
    names_to = "model",
    values_to = "return_level",
    names_prefix = "levels_"
  ) |>
  mutate(
    model = case_when(
      model == "forest" ~ "Random Forest",
      model == "gp" ~ "GP conversion",
      TRUE ~ model
    )
  )

write_rds(
  return_levels_long,
  here::here("data", "era5_land_hourly_alps_dl_return_levels.rds")
)

log_info("Finished 5-runoff_marginals.r")
