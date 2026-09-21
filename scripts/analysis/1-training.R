# Stage 1: build the training table (predictors at POT peak hours) and the
# T-year runoff event levels.
# Reads:  derived/era5_land_hourly_alps_{all,peaks}.rds (scripts 2-3)
# Writes: out/training.rds, out/event_levels.rds
source(here::here("scripts", "lib", "common.R"))

peaks <- read_rds(here::here("derived", "era5_land_hourly_alps_peaks.rds"))
hourly <- read_rds(here::here("derived", "era5_land_hourly_alps_all.rds"))

log_info("Computing predictors: {paste(cfg$predictors, collapse = ', ')}")
training <- build_training(hourly, peaks, cfg$predictors, cfg$response)
n_missing <- sum(!complete.cases(training))
if (n_missing > 0) log_warn("{n_missing} peak rows have missing predictors")
write_rds(training, out_path("training.rds"))

log_info("T-year event levels from a GP fit to each cell's POT peaks")
n_years <- diff(range(year(hourly$date))) + 1
rm(hourly)
rp_curve <- sort(unique(c(cfg$queries$return_periods, exp(seq(log(1.1), log(500), length.out = 80)))))
event_levels <- peaks |>
  group_by(cell_id, x, y) |>
  group_modify(\(d, key) {
    pot_event_levels(d$runoff_hourly, d$threshold[1], n_years, rp_curve)
  }) |>
  ungroup() |>
  mutate(in_queries = return_period %in% cfg$queries$return_periods)
write_rds(event_levels, out_path("event_levels.rds"))
log_info("Stage 1 done: {nrow(training)} peaks over {n_distinct(training$cell_id)} cells")
