# Plotting helpers shared by scripts/analysis/5-figures.R and apps/explorer.
# Every function takes stage outputs (tibbles / models) and returns a ggplot.

VAR_LABELS <- c(
  runoff_hourly = "Runoff (mm/h)",
  rainfall_hourly = "Rainfall (mm/h)",
  snowmelt_hourly = "Snowmelt (mm/h)",
  swe_mm = "Snow water equivalent (mm)",
  rainfall_24h = "Rainfall, past 24 h (mm)",
  snowmelt_24h = "Snowmelt, past 24 h (mm)"
)
var_label <- function(v) ifelse(v %in% names(VAR_LABELS), VAR_LABELS[v], v)

# Ordinal blue ramp (light -> dark) for return periods; grey for context.
BLUE_RAMP <- c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b")
ORANGE <- "#eb6834"
INK <- "#0b0b0b"
INK_2 <- "#52514e"
GRID <- "#e6e5e0"

rp_colours <- function(rps) {
  rps <- sort(unique(rps))
  cols <- grDevices::colorRampPalette(BLUE_RAMP)(max(2, length(rps)))
  stats::setNames(cols[seq_along(rps)], paste0(rps, "-yr"))
}
rp_factor <- function(rp) factor(paste0(rp, "-yr"), levels = paste0(sort(unique(rp)), "-yr"))

theme_ros <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = GRID, linewidth = 0.3),
      axis.title = ggplot2::element_text(colour = INK_2),
      axis.text = ggplot2::element_text(colour = INK_2),
      plot.title = ggplot2::element_text(colour = INK, face = "bold"),
      plot.subtitle = ggplot2::element_text(colour = INK_2),
      strip.text = ggplot2::element_text(colour = INK, face = "bold", hjust = 0),
      legend.title = ggplot2::element_text(colour = INK_2),
      legend.position = "top",
      legend.justification = "left",
      plot.title.position = "plot"
    )
}

cell_label <- function(cell_id, x, y) sprintf("Cell %s (%.1f°E, %.1f°N)", cell_id, x, y)

# ---- Data ----------------------------------------------------------------

plot_training_scatter <- function(training, cfg) {
  p <- cfg$predictors
  d <- dplyr::mutate(training, cell = cell_label(cell_id, x, y))
  if (length(p) == 1) {
    ggplot2::ggplot(d, ggplot2::aes(.data[[p]], .data[[cfg$response]])) +
      ggplot2::geom_point(alpha = 0.35, size = 1.2, colour = INK_2) +
      ggplot2::facet_wrap(~cell) +
      ggplot2::labs(x = var_label(p), y = var_label(cfg$response),
        title = "Runoff peaks against rainfall", subtitle = "One point per POT peak hour") +
      theme_ros()
  } else {
    ggplot2::ggplot(d, ggplot2::aes(.data[[p[1]]], .data[[p[2]]], colour = .data[[cfg$response]])) +
      ggplot2::geom_point(alpha = 0.7, size = 1.2) +
      ggplot2::scale_colour_gradient(low = "#cde2fb", high = "#0d366b", name = var_label(cfg$response)) +
      ggplot2::facet_wrap(~cell) +
      ggplot2::labs(x = var_label(p[1]), y = var_label(p[2]),
        title = "Drivers at runoff peaks", subtitle = "Colour = runoff at the peak hour") +
      theme_ros()
  }
}

# ---- LOESS-like model, one predictor ---------------------------------------

# Scatter with fitted conditional quantile curves across the predictor range.
plot_llqr_curves <- function(model, training, levels = c(0.1, 0.5, 0.9, 0.99), n = 80,
                             title = "Fitted conditional quantiles") {
  xn <- model$xnames; yn <- model$yname
  xs <- seq(min(training[[xn]]), max(training[[xn]]), length.out = n)
  q <- llqr_quantiles(model, xs, levels = levels)
  curves <- tibble::tibble(x0 = rep(xs, times = length(levels)),
    level = rep(levels, each = n), q = as.numeric(q))
  labs_end <- dplyr::filter(curves, x0 == max(x0))
  cols <- grDevices::colorRampPalette(BLUE_RAMP[c(1, 3, 5, 6)])(length(levels))
  ggplot2::ggplot() +
    ggplot2::geom_point(data = training, ggplot2::aes(.data[[xn]], .data[[yn]]),
      colour = INK_2, alpha = 0.25, size = 1.1) +
    ggplot2::geom_line(data = curves, ggplot2::aes(x0, q, colour = factor(level)), linewidth = 1) +
    ggplot2::geom_text(data = labs_end, ggplot2::aes(x0, q, label = paste0("τ = ", level)),
      hjust = -0.1, size = 3.3, colour = INK_2) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.15))) +
    ggplot2::labs(x = var_label(xn), y = var_label(yn), title = title,
      subtitle = sprintf("Local linear quantile regression, span = %s, %s kernel", model$span, model$kernel)) +
    theme_ros()
}

# How one prediction is made: weights around x0, local lines, read-off at x0.
plot_llqr_local <- function(model, x0, levels = c(0.1, 0.5, 0.9)) {
  xn <- model$xnames; yn <- model$yname
  tr <- model$training
  x <- tr[[xn]]; y <- tr[[yn]]
  h <- model$bandwidth %||% rainonsnow:::llqr_nn_bandwidth(x, x0, model$span)
  w <- llqr_weights(x, x0, h, model$kernel)
  lines <- purrr::map_dfr(levels, function(p) {
    keep <- w > 0
    ab <- llqr_local_line(x[keep] - x0, y[keep], w[keep], p, model$degree)
    xs <- seq(max(min(x), x0 - h), min(max(x), x0 + h), length.out = 30)
    tibble::tibble(level = p, xx = xs, yy = ab[1] + ab[2] * (xs - x0), a = ab[1])
  })
  at_x0 <- dplyr::distinct(lines, level, a)
  pts <- dplyr::mutate(tr, weight = w)
  ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = x0 - h, xmax = x0 + h, ymin = -Inf, ymax = Inf,
      fill = "#cde2fb", alpha = 0.35) +
    ggplot2::geom_point(data = dplyr::filter(pts, weight == 0),
      ggplot2::aes(.data[[xn]], .data[[yn]]), colour = "#c3c2b7", size = 0.9, alpha = 0.5) +
    ggplot2::geom_point(data = dplyr::filter(pts, weight > 0),
      ggplot2::aes(.data[[xn]], .data[[yn]], alpha = weight, size = weight), colour = "#1c5cab") +
    ggplot2::geom_vline(xintercept = x0, colour = INK_2, linetype = "dashed") +
    ggplot2::geom_line(data = lines, ggplot2::aes(xx, yy, group = level), colour = ORANGE, linewidth = 1) +
    ggplot2::geom_point(data = at_x0, ggplot2::aes(x0, a), colour = ORANGE, size = 3.2) +
    ggplot2::geom_label(data = at_x0, ggplot2::aes(x0, a, label = paste0("τ=", level)),
      hjust = 1.15, size = 3, label.size = 0, fill = "white", colour = INK) +
    ggplot2::scale_alpha(range = c(0.15, 0.9), guide = "none") +
    ggplot2::scale_size(range = c(0.6, 2.6), guide = "none") +
    ggplot2::labs(x = var_label(xn), y = var_label(yn),
      title = sprintf("Local fit at %s = %.2f", var_label(xn), x0),
      subtitle = "Shaded window: observations weighted by distance from x0 (darker = more weight).\nOrange: weighted quantile lines; the dot is the conditional quantile at x0.") +
    theme_ros()
}

# Single-level local line (intercept, slope) at x0, for illustration.
llqr_local_line <- function(xc, y, w, p, degree) {
  a <- rainonsnow:::llqr_solve_one(xc, y, w, p, degree)
  c(a, attr(a, "slope"))
}

# ---- Predictive distributions -------------------------------------------

# Exceedance curves P(runoff > z | x) at several predictor settings (rows of
# newdata), with the T-year event levels marked.
plot_predictive_exceedance <- function(model, newdata, tail, event_levels) {
  dsts <- dl_apply_tail(predict(model, newdata), tail)
  zmax <- max(event_levels$return_level, na.rm = TRUE) * 1.15
  zs <- seq(0, zmax, length.out = 300)
  lab <- apply(newdata[model$xnames], 1, function(r) {
    paste(sprintf("%s = %.2f", sub("_hourly", "", model$xnames), r), collapse = ", ")
  })
  curves <- purrr::map2_dfr(dsts, seq_along(dsts), function(d, i) {
    if (identical(distionary::pretty_name(d), "Null")) return(NULL)
    tibble::tibble(i = i, setting = lab[i], z = zs,
      surv = as.numeric(distionary::eval_survival(d, at = zs)))
  })
  curves$setting <- factor(curves$setting, levels = unique(lab))
  lv <- dplyr::filter(event_levels, in_queries)
  cols <- grDevices::colorRampPalette(BLUE_RAMP)(max(2, nrow(newdata)))[seq_len(nrow(newdata))]
  ggplot2::ggplot(dplyr::filter(curves, surv > 1e-4), ggplot2::aes(z, surv, colour = setting)) +
    ggplot2::geom_vline(data = lv, ggplot2::aes(xintercept = return_level), colour = GRID, linewidth = 0.6) +
    ggplot2::geom_text(data = lv, ggplot2::aes(x = return_level, y = 1, label = paste0(return_period, "-yr")),
      inherit.aes = FALSE, angle = 90, hjust = 1, vjust = -0.4, size = 3, colour = INK_2) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_y_log10(labels = scales::label_percent(drop0trailing = TRUE)) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
    ggplot2::labs(x = var_label(model$yname), y = "P(runoff > z | drivers)",
      title = "Predictive distributions, with GP tails",
      subtitle = "Chance of exceeding each runoff level; vertical lines mark the T-year levels") +
    theme_ros()
}

# ---- Queries -------------------------------------------------------------

plot_event_curve_1d <- function(event, triggers, target, prob = 0.5) {
  ev <- dplyr::mutate(event, T = rp_factor(return_period))
  tg <- dplyr::filter(triggers, abs(.data$prob - .env$prob) < 1e-9, !is.na(trigger)) |>
    dplyr::mutate(T = rp_factor(return_period))
  cols <- rp_colours(event$return_period)
  ggplot2::ggplot(ev, ggplot2::aes(.data[[target]], p_event, colour = T)) +
    ggplot2::geom_hline(yintercept = prob, colour = INK_2, linetype = "dotted") +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(data = tg, ggplot2::aes(trigger, .env$prob), size = 3, shape = 21,
      fill = "white", stroke = 1.3) +
    ggplot2::scale_colour_manual(values = cols, name = "Event") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = var_label(target), y = "P(runoff exceeds the T-year level)",
      title = "Chance that rainfall triggers a T-year runoff peak",
      subtitle = sprintf("Open circles: rainfall at which the chance reaches %s%%", prob * 100)) +
    theme_ros()
}

plot_trigger_2d <- function(triggers, target, other, rps = NULL) {
  tg <- triggers
  if (!is.null(rps)) tg <- dplyr::filter(tg, return_period %in% rps)
  tg <- dplyr::mutate(tg,
    T = rp_factor(return_period),
    chance = factor(paste0(round(prob * 100), "% chance"), levels = paste0(sort(unique(round(prob * 100))), "% chance"))
  )
  ggplot2::ggplot(tg, ggplot2::aes(.data[[other]], trigger, colour = chance, linewidth = chance)) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::facet_wrap(~T, labeller = ggplot2::labeller(T = function(x) paste(x, "event"))) +
    ggplot2::scale_colour_manual(values = c("#86b6ef", "#1c5cab", "#0d366b"), name = NULL) +
    ggplot2::scale_linewidth_manual(values = c(0.7, 1.4, 0.7), name = NULL) +
    ggplot2::labs(x = var_label(other), y = paste("Rain needed:", var_label(target)),
      title = "Rain needed for a T-year runoff peak, given snowmelt",
      subtitle = "Gaps: that chance is not reached anywhere in the observed rainfall range") +
    theme_ros()
}

plot_event_surface_2d <- function(event, training, target, other, rp) {
  ev <- dplyr::filter(event, return_period == rp)
  ggplot2::ggplot(ev, ggplot2::aes(.data[[target]], .data[[other]])) +
    ggplot2::geom_raster(ggplot2::aes(fill = p_event), interpolate = TRUE) +
    ggplot2::geom_contour(ggplot2::aes(z = p_event), breaks = c(0.1, 0.5, 0.9), colour = "white", linewidth = 0.5) +
    ggplot2::geom_point(data = training, colour = INK, alpha = 0.3, size = 0.6) +
    ggplot2::scale_fill_gradientn(colours = c("#f6f7f5", "#cde2fb", "#86b6ef", "#2a78d6", "#0d366b"),
      labels = scales::percent, name = "P(event)", limits = c(0, 1)) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::labs(x = var_label(target), y = var_label(other),
      title = sprintf("Chance of a %s-year runoff peak", rp),
      subtitle = "White contours: 10%, 50%, 90%. Dots: observed peaks.") +
    theme_ros() + ggplot2::theme(legend.position = "right")
}

plot_likeliest_1d <- function(likeliest, target) {
  lk <- dplyr::mutate(likeliest, T = rp_factor(return_period))
  prior <- dplyr::distinct(likeliest, .data[[target]], f_prior)
  cols <- rp_colours(likeliest$return_period)
  ggplot2::ggplot(lk, ggplot2::aes(.data[[target]], f_given_event, colour = T)) +
    ggplot2::geom_area(data = prior, ggplot2::aes(.data[[target]], f_prior), inherit.aes = FALSE,
      fill = "#e6e5e0", colour = NA) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_manual(values = cols, name = "Given a") +
    ggplot2::labs(x = var_label(target), y = "Density",
      title = "What rainfall is behind a T-year runoff peak?",
      subtitle = "Grey: rainfall at all runoff peaks. Lines: rainfall given the T-year level was exceeded.") +
    theme_ros()
}

plot_likeliest_2d_cond <- function(likeliest, target, other, rp, n_slices = 4) {
  lk <- dplyr::filter(likeliest, return_period == rp)
  vals <- sort(unique(lk[[other]]))
  pick <- vals[unique(round(seq(1, length(vals) * 0.8, length.out = n_slices)))]
  lk <- dplyr::filter(lk, .data[[other]] %in% pick) |>
    dplyr::mutate(slice = factor(sprintf("%.2f", .data[[other]])))
  cols <- grDevices::colorRampPalette(BLUE_RAMP[c(1, 3, 5, 6)])(length(pick))
  ggplot2::ggplot(lk, ggplot2::aes(.data[[target]], f_given_event_cond, colour = slice)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_colour_manual(values = cols, name = var_label(other)) +
    ggplot2::labs(x = var_label(target), y = "Density",
      title = sprintf("Rainfall behind a %s-year runoff peak, by snowmelt", rp),
      subtitle = "Density of rainfall given the event and the snowmelt rate") +
    theme_ros()
}

# ---- Diagnostics ---------------------------------------------------------

plot_return_levels <- function(return_levels) {
  d <- dplyr::mutate(return_levels, cell = cell_label(cell_id, x, y))
  ggplot2::ggplot(d, ggplot2::aes(return_period, return_level, colour = source, linetype = source)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200, 500)) +
    ggplot2::scale_colour_manual(values = c("#2a78d6", ORANGE), name = NULL) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed"), name = NULL) +
    ggplot2::facet_wrap(~cell, scales = "free_y") +
    ggplot2::labs(x = "Return period (years)", y = "Hourly runoff peak (mm/h)",
      title = "Runoff return levels",
      subtitle = "Does the DL model reproduce the peaks it was trained on?") +
    theme_ros()
}

plot_pp <- function(diagnostics) {
  d <- dplyr::mutate(diagnostics$pp, model = dplyr::recode(model, raw = "Raw prediction", gp = "With GP tail"))
  ggplot2::ggplot(d, ggplot2::aes(p_empirical, p_model, group = cell_id)) +
    ggplot2::geom_abline(colour = INK_2, linetype = "dashed") +
    ggplot2::geom_line(colour = "#2a78d6", alpha = 0.7) +
    ggplot2::facet_wrap(~model) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Empirical probability", y = "Model PIT",
      title = "Calibration (in-sample P-P)", subtitle = "One line per cell; on the diagonal = calibrated") +
    theme_ros()
}

plot_skill <- function(diagnostics) {
  ggplot2::ggplot(diagnostics$skill, ggplot2::aes(tau, skill_score, group = cell_id)) +
    ggplot2::geom_hline(yintercept = 0, colour = INK_2) +
    ggplot2::geom_line(colour = "#2a78d6", alpha = 0.8) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(x = "Quantile level", y = "Skill vs. unconditional",
      title = "Quantile skill score", subtitle = "Improvement in pinball loss over the cell's runoff marginal (in-sample)") +
    theme_ros()
}
