# Stage 4: query the fitted joint distribution of (runoff, predictors).
#   - event probability  P(runoff > z_T | x) over a predictor grid
#   - trigger thresholds  smallest target predictor (e.g. rain) giving
#                         P(event) >= prob, for each value of the others
#   - likeliest drivers  f(x | runoff > z_T), and with two predictors
#                         f(target | runoff > z_T, other)
# Writes: out/drivers.rds, out/event_probability.rds, out/triggers.rds,
#         out/likeliest.rds
source(here::here("scripts", "lib", "common.R"))

training <- read_rds(out_path("training.rds"))
models <- read_rds(out_path("models.rds"))
event_levels <- read_rds(out_path("event_levels.rds")) |> filter(in_queries)
target <- cfg$queries$trigger$target
log_info("Queries for target '{target}', T = {paste(cfg$queries$return_periods, collapse = ', ')} years")

res <- map_cells(training, \(d) {
  d <- tidyr::drop_na(d)
  id <- d$cell_id[1]
  model <- models$model[[match(id, models$cell_id)]]
  lv <- filter(event_levels, cell_id == id)
  grid <- query_grid(d, cfg$predictors, cfg$queries$grid$n, cfg$queries$grid$mult)
  ev <- query_event_probability(model, grid, lv, cfg$tail)
  drivers <- fit_drivers(d, cfg$predictors, cfg$drivers$marginal_family, cfg$drivers$family_set)
  keys <- tibble(cell_id = id, x = d$x[1], y = d$y[1])
  list(
    drivers = mutate(keys, drivers = list(drivers)),
    event = bind_cols(keys, ev),
    trigger = bind_cols(keys, query_trigger(ev, target, cfg$queries$trigger$probs)),
    likeliest = bind_cols(keys, query_likeliest(ev, drivers, target))
  )
})

write_rds(bind_rows(map(res, "drivers")), out_path("drivers.rds"))
write_rds(bind_rows(map(res, "event")), out_path("event_probability.rds"))
write_rds(bind_rows(map(res, "trigger")), out_path("triggers.rds"))
write_rds(bind_rows(map(res, "likeliest")), out_path("likeliest.rds"))
log_info("Stage 4 done")
