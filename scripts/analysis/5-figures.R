# Stage 5: figures for the analysis README, the app and slides.
# Writes: figs/*.png
source(here::here("scripts", "lib", "common.R"))
source(here::here("scripts", "lib", "plots.R"))

training <- read_rds(out_path("training.rds"))
models <- read_rds(out_path("models.rds"))
event_levels <- read_rds(out_path("event_levels.rds"))
diagnostics <- read_rds(out_path("diagnostics.rds"))
return_levels <- read_rds(out_path("return_levels.rds"))
event <- read_rds(out_path("event_probability.rds"))
triggers <- read_rds(out_path("triggers.rds"))
likeliest <- read_rds(out_path("likeliest.rds"))

focus <- choose_focus_cell(cfg, training)
log_info("Focus cell: {focus}")
fc <- function(d) filter(d, cell_id == focus)
model <- models$model[[match(focus, models$cell_id)]]
tr_focus <- fc(training)
focus_lab <- with(tr_focus[1, ], cell_label(cell_id, x, y))
target <- cfg$queries$trigger$target

save_fig <- function(p, name, w = 8, h = 5.5) {
  ggsave(fig_path(paste0(name, ".png")), p, width = w, height = h, dpi = 200, bg = "white")
  log_info("Wrote figs/{name}.png")
}
add_cell <- function(p) p + labs(caption = focus_lab)

save_fig(plot_training_scatter(training, cfg), "01-training-scatter", 9, 7)
save_fig(plot_return_levels(return_levels), "02-return-levels", 9, 7)
save_fig(plot_pp(diagnostics), "03-calibration-pp", 9, 5)
save_fig(plot_skill(diagnostics), "04-quantile-skill", 8, 5)

if (length(cfg$predictors) == 1) {
  xq <- quantile(tr_focus[[target]], c(0.25, 0.75, 0.97))
  if (cfg$model$type == "llqr") {
    save_fig(add_cell(plot_llqr_curves(model, tr_focus)), "10-llqr-quantile-curves")
    for (i in seq_along(xq)) {
      save_fig(add_cell(plot_llqr_local(model, unname(xq[i]))), sprintf("11-llqr-local-fit-%d", i))
    }
  }
  xs <- signif(seq(0.5, max(tr_focus[[target]]), length.out = 6), 2)
  nd <- tibble(!!target := xs)
  save_fig(add_cell(plot_predictive_exceedance(model, nd, cfg$tail, fc(event_levels))),
    "12-predictive-exceedance")
  save_fig(add_cell(plot_event_curve_1d(fc(event), fc(triggers), target)), "20-event-probability")
  save_fig(add_cell(plot_likeliest_1d(fc(likeliest), target)), "21-likeliest-rainfall")
} else {
  other <- setdiff(cfg$predictors, target)
  save_fig(add_cell(plot_trigger_2d(fc(triggers), target, other, rps = c(2, 5, 10))),
    "20-trigger-rain-given-snowmelt", 10, 5.5)
  for (rp in c(2, 10, 50)) {
    save_fig(add_cell(plot_event_surface_2d(fc(event), tr_focus, target, other, rp)),
      sprintf("21-event-surface-T%03d", rp), 8, 6)
  }
  save_fig(add_cell(plot_likeliest_2d_cond(fc(likeliest), target, other, 10)),
    "22-likeliest-rain-given-snowmelt")
  save_fig(add_cell(plot_likeliest_2d_cond(fc(likeliest), target, other, 2)),
    "22-likeliest-rain-given-snowmelt-T002")
}
log_info("Stage 5 done")
