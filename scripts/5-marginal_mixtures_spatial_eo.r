# Marginal distribution modelling: equal-weight mixture of peak-hour predictive
# distributions per cell (forest + GP tail) and precomputed return-level table
# for apps/return-level-explorer.
# Requires: outputs of scripts/3-pot_spatial_eo.r and 4-dl_rqforest_spatial_eo.r
# %%
library(tidyverse)
library(logger)
devtools::load_all()

return_periods <- rp_reporting()

log_info("Starting 5-marginal_mixtures_spatial_eo.r")

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

lvls <- marginals_tbl |>
  mutate(
    df = map(
      marginal_forest,
      enframe_return,
      at = return_periods,
      arg_name = "return_period"
    )
  ) |>
  select(cell_id, x, y, df) |>
  unnest(df)

lvlsgp <- marginals_tbl |>
  left_join(num_pot_events, by = c("cell_id", "x", "y")) |>
  mutate(
    df = map2(
      marginal_gp,
      num_events_per_year,
      \(mr, nep) {
        enframe_return(mr, at = return_periods * nep, arg_name = "rp_adjusted")
      }
    )
  ) |>
  select(cell_id, x, y, df) |>
  unnest(df) |>
  mutate(return_period = rp_adjusted / num_events_per_year)

lvls_both <- left_join(
  lvls,
  lvlsgp |> select(cell_id, x, y, return_period, return),
  by = c("cell_id", "x", "y", "return_period"),
  suffix = c("_emp", "_gp")
)

lvls_both_long <- lvls_both |>
  pivot_longer(
    c(return_emp, return_gp),
    names_to = "model",
    values_to = "return_level",
    names_prefix = "return_"
  ) |>
  mutate(
    model = case_when(
      model == "emp" ~ "Forest mixture (empirical)",
      model == "gp" ~ "Forest mixture + GP tail",
      TRUE ~ model
    )
  )

write_rds(
  list(lvls_both_long = lvls_both_long, return_periods = return_periods),
  here::here("data", "era5_land_hourly_alps_dl_marginal_return_levels.rds")
)

log_info("Finished 5-marginal_mixtures_spatial_eo.r")
