# Stage 2: fit the distributional learning model per cell, predict at every
# peak (raw and tail-treated), and compute calibration diagnostics.
# Writes: out/models.rds, out/predictions.rds (compact), out/diagnostics.rds
source(here::here("scripts", "lib", "common.R"))

training <- read_rds(out_path("training.rds"))
log_info("Fitting '{cfg$model$type}' on {N_CORES} core(s)")

fits <- map_cells(training, \(d) {
  d <- tidyr::drop_na(d)
  model <- dl_fit(d, cfg)
  raw <- predict(model, newdata = d)
  list(
    model = model,
    predictions = mutate(d,
      distribution_raw = raw,
      distribution_gp = dl_apply_tail(raw, cfg$tail)
    )
  )
})

models <- training |>
  distinct(cell_id, x, y) |>
  arrange(cell_id) |>
  mutate(model = map(as.character(cell_id), \(id) fits[[id]]$model))
write_rds(models, out_path("models.rds"))

predictions <- bind_rows(map(fits, "predictions"))
n_null <- sum(map_chr(predictions$distribution_gp, pretty_name) == "Null")
if (n_null > 0) log_warn("{n_null} peak predictions have no usable tail")

dl_write_peak_hour_predictions(predictions, out_path("predictions.rds"))

log_info("Diagnostics (P-P calibration, quantile skill vs marginal)")
diag_input <- predictions |>
  filter(map_chr(distribution_gp, pretty_name) != "Null")
diagnostics <- dl_build_diagnostics(diag_input, training)
write_rds(diagnostics, out_path("diagnostics.rds"))
log_info("Stage 2 done")
