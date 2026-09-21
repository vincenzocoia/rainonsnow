# Stage 3: runoff return levels implied by the DL model (mixture of peak-hour
# predictive distributions), alongside the direct POT fit from stage 1. The
# comparison checks that the DL model reproduces the marginal it was trained on.
# Writes: out/return_levels.rds
source(here::here("scripts", "lib", "common.R"))

predictions <- dl_read_peak_hour_predictions(out_path("predictions.rds"))
event_levels <- read_rds(out_path("event_levels.rds"))

dl_levels <- map_cells(predictions, \(d) {
  rate <- event_levels$rate[event_levels$cell_id == d$cell_id[1]][1]
  rps <- unique(event_levels$return_period)
  rps <- rps[rps * rate > 1]
  tibble(
    cell_id = d$cell_id[1], x = d$x[1], y = d$y[1],
    return_period = rps,
    return_level = mixture_return_levels(d$distribution_gp, rps, rate),
    source = "DL model (mixture over peaks)"
  )
})

return_levels <- bind_rows(
  event_levels |>
    transmute(cell_id, x, y, return_period, return_level, source = "POT peaks (GP fit)"),
  bind_rows(dl_levels)
)
write_rds(return_levels, out_path("return_levels.rds"))
log_info("Stage 3 done")
