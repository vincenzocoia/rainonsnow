# Same conditional surface as 7-likeliest_rain_snow-one_cell_animated.r, but with
# BOTH driver axes on their own return-period scale.
#
# THE QUESTION THIS ANSWERS
#
# The raw-axis version says the likeliest driver pair for a T-year runoff event
# is (6.8 mm/h rain, 0.3 mm/h snowmelt). That is not readable as a statement
# about rarity: is 0.3 mm/h of snowmelt a lot? Here each axis is relabelled as
# "an event of this size happens once every T_R (or T_S) years at this cell", so
# a point on the surface reads directly as "a T_R-year rainfall together with a
# T_S-year snowmelt produces the T-year runoff".
#
# WHY THE DRIVER MARGINALS ARE FITTED HERE AND NOT REUSED
#
# derived/..._peaks.rds holds rainfall and snowmelt AT RUNOFF PEAKS. Those are
# conditional on runoff being extreme, so their empirical distribution is not the
# marginal distribution of rainfall or of snowmelt, and a return period read off
# them would be meaningless. Each driver therefore gets its own peaks-over-
# threshold analysis on the FULL hourly series for the cell, declustered on its
# own timing, with the same settings script 3 uses for runoff
# (inputs/pot_metadata.yaml).
#
# WHY THE TAIL ALONE IS NOT ENOUGH
#
# A GP tail above the 99th percentile only defines return periods above roughly
# 1/lambda years. For rainfall that is fine -- the interesting region is far out
# in the tail. For snowmelt it is not: the snowmelt accompanying even a 200-year
# runoff event sits BELOW its own POT threshold, so a pure GPD would leave the
# whole surface off the bottom of the axis. Each driver marginal is therefore
# empirical in the body and generalized Pareto in the tail, which is the same
# body/tail split script 5 already uses for the runoff marginal.
#
# THE SURFACE IS TRANSFORMED, NOT JUST RELABELLED
#
# T_R is a monotone function of rainfall, so moving to return-period axes is a
# change of variable, and a density has to carry the Jacobian with it. What is
# drawn is the density of (log T_R, log T_S) given the runoff level, so areas on
# the plot are probability. The mode therefore need not sit at the transform of
# the raw-axis mode; both are reported.
#
# Requires everything script 7 requires, plus:
#   - derived/era5_land_hourly_alps_all.rds (script 2) for the driver POT fits
# %%
library(tidyverse)
library(rvinecopulib)
library(yaml)
library(magick)
devtools::load_all()

# --- edit here ---
# Cell to plot: a command-line argument wins, then the config, then auto-select.
# The argument is what lets this run over every cell:
#   for c in 1 2 3 4; do Rscript <this script> $c; done
.arg <- commandArgs(trailingOnly = TRUE)
.cfg <- read_yaml(here::here("inputs", "rain_snow_joint_model.yaml"))$likeliest_rain_snow
cell_id <- suppressWarnings(as.integer(
  if (length(.arg) >= 1L) .arg[[1]] else .cfg$cell_id %||% NA_integer_
))
if (is.na(cell_id)) {
  .peaks_all <- read_rds(here::here("derived", "era5_land_hourly_alps_peaks.rds"))
  .mixed <- .peaks_all |>
    filter(rainfall_hourly > 0, snowmelt_hourly > 0) |>
    count(cell_id, sort = TRUE)
  cell_id <- if (nrow(.mixed) > 0) {
    .mixed$cell_id[[1]]
  } else {
    sort(unique(.peaks_all$cell_id))[[1]]
  }
  message("Focus cell auto-selected (most mixed rain+snow peaks): ", cell_id)
}

grid_size <- c(45L, 45L)
still_targets <- c(2, 10, 50, 200)
return_period_years <- exp(seq(log(2), log(200), length.out = 30))
for (.t in still_targets) {
  return_period_years[which.min(abs(log(return_period_years) - log(.t)))] <- .t
}
return_period_years <- sort(return_period_years)
stopifnot(!anyDuplicated(sprintf("%g", return_period_years)))
fps <- 4

# Driver POT settings: the same ones script 3 applies to runoff.
.pot <- read_yaml(here::here("inputs", "pot_metadata.yaml"))
pot_quantile <- as.numeric(.pot$quantile %||% 0.99)
pot_min_gap <- as.integer(.pot$min_gap %||% 72)
# Declustering threshold for the BODY of each driver marginal. Lower than the
# POT threshold so the marginal reaches return periods below 1/lambda years,
# which is where the snowmelt axis has to reach.
body_quantile <- 0.90

out_gif <- here::here(
  "plots",
  sprintf("rain_snow_return_period_axes_cell_%d_animated.gif", cell_id)
)
driver_cache <- here::here(
  "derived",
  "era5_land_hourly_alps_driver_pot_marginals.rds"
)
# ---

# %%
peaks_raw <- read_rds(here::here("derived", "era5_land_hourly_alps_peaks.rds"))
peaks <- filter(peaks_raw, rainfall_hourly != 0, snowmelt_hourly != 0)
peaks_cell <- filter(peaks, cell_id == .env$cell_id)

# The largest runoff this cell has actually produced in the record. Any frame
# whose return level exceeds it is extrapolation, and a frame far above it is a
# warning that the cell's mixture tail is being driven by one heavy component
# rather than by its data -- see the check below.
max_obs_runoff <- max(filter(peaks_raw, cell_id == .env$cell_id)$runoff_hourly)

joint <- read_rds(
  here::here("derived", "era5_land_hourly_alps_joint_rain_snow.rds")
) |>
  filter(cell_id == .env$cell_id) |>
  pull(joint) |>
  pluck(1)

dl_model <- read_rds(
  here::here("derived", "era5_land_hourly_alps_dl_rqforest_models.rds")
) |>
  filter(cell_id == .env$cell_id) |>
  pull(dl_rqforest) |>
  pluck(1)

tail_row <- read_rds(
  here::here("derived", "era5_land_hourly_alps_dl_mixture_tails.rds")
) |>
  filter(cell_id == .env$cell_id)
cell_mt <- tail_row$mixture_tail[[1]]
nep <- tail_row$num_events_per_year[1]

cell_shape <- read_rds(
  here::here("derived", "era5_land_hourly_alps_dl_tail_shapes.rds")
) |>
  filter(cell_id == .env$cell_id) |>
  pull(shape_raw)
stopifnot(length(cell_shape) == 1L, is.finite(cell_shape))

# %%
# ---- Driver marginals: an independent POT analysis per driver ----------------

#' One driver's marginal, empirical in the body and generalized Pareto above the
#' POT threshold, expressed as a return period in years.
fit_driver_marginal <- function(df, col, n_years) {
  u_body <- unname(stats::quantile(df[[col]], body_quantile, na.rm = TRUE))
  events <- get_pot_events(
    df,
    threshold = u_body,
    flow_col = col,
    date_col = "date",
    min_gap = pot_min_gap
  )
  x <- sort(events[[col]])
  n <- length(x)
  lambda <- n / n_years

  u <- unname(stats::quantile(df[[col]], pot_quantile, na.rm = TRUE))
  excess <- x[x > u] - u
  shape_bounds <- c(-0.45, 1)
  fit <- fit_gpd_weighted(excess, rep(1, length(excess)), shape_bounds = shape_bounds)
  zeta_u <- mean(x > u)

  # Weibull plotting positions give a smooth monotone body; the GP takes over at
  # u, where the two agree by construction (the GP survival is 1 at u).
  p_emp <- (n - seq_len(n) + 1) / (n + 1)

  surv <- function(v) {
    out <- numeric(length(v))
    lo <- v <= u
    out[lo] <- stats::approx(x, p_emp, xout = v[lo], rule = 2)$y
    out[!lo] <- zeta_u *
      gpd_survival(v[!lo] - u, fit$scale, fit$shape)
    pmax(out, .Machine$double.eps)
  }

  list(
    col = col,
    n_events = n,
    n_years = n_years,
    lambda = lambda,
    u_body = u_body,
    u_pot = u,
    n_exceed = length(excess),
    gp_scale = fit$scale,
    gp_shape = fit$shape,
    zeta_u = zeta_u,
    x_min = min(x),
    x_max = max(x),
    events = x,
    # fit_gpd_weighted() searches xi over a bounded interval. A shape sitting ON
    # a bound is not an optimum, it is the search giving up, and everything the
    # tail says above the threshold is then unsupported.
    # stats::optimize() never returns an endpoint exactly; its default tolerance
    # is .Machine$double.eps^0.25 (~1.2e-4), so "on the bound" has to be judged
    # at that scale, not at 1e-6.
    shape_at_bound = min(abs(fit$shape - shape_bounds)) <
      10 * .Machine$double.eps^0.25,
    # T in years for a level v: one event above v every 1 / (lambda * S(v)) years.
    T_of_x = function(v) 1 / (lambda * surv(v))
  )
}

if (fs::file_exists(driver_cache)) {
  driver_fits <- read_rds(driver_cache)
} else {
  driver_fits <- list()
}
cache_key <- as.character(cell_id)

# The hourly file is ~300 MB, so a cache miss fits EVERY cell while it is open
# rather than paying that read once per cell.
all_cells <- sort(unique(read_rds(
  here::here("derived", "era5_land_hourly_alps_peaks.rds")
)$cell_id))
missing_cells <- setdiff(as.character(all_cells), names(driver_fits))

if (length(missing_cells) > 0L) {
  message(
    "Reading the full hourly series to fit driver marginals for cell(s) ",
    paste(missing_cells, collapse = ", "),
    " ..."
  )
  hourly_all <- read_rds(here::here("derived", "era5_land_hourly_alps_all.rds")) |>
    select(cell_id, date, rainfall_hourly, snowmelt_hourly)

  for (ck in missing_cells) {
    h <- filter(hourly_all, cell_id == as.integer(ck))
    n_years <- diff(range(year(h$date))) + 1
    driver_fits[[ck]] <- list(
      rain = fit_driver_marginal(h, "rainfall_hourly", n_years),
      snow = fit_driver_marginal(h, "snowmelt_hourly", n_years)
    )
    message("  fitted cell ", ck)
  }
  rm(hourly_all)
  gc(verbose = FALSE)
  write_rds(driver_fits, driver_cache)
  message("Cached driver marginals to ", driver_cache)
} else {
  message("Using cached driver marginals for cell ", cell_id, ".")
}
stopifnot(!is.null(driver_fits[[cache_key]]))

rain_m <- driver_fits[[cache_key]]$rain
snow_m <- driver_fits[[cache_key]]$snow

for (m in list(rain_m, snow_m)) {
  message(sprintf(
    paste0(
      "%s: %d declustered events over %d y (%.2f/yr); POT threshold %.3f mm/h ",
      "with %d exceedances, GP scale %.3f shape %.3f."
    ),
    m$col,
    m$n_events,
    m$n_years,
    m$lambda,
    m$u_pot,
    m$n_exceed,
    m$gp_scale,
    m$gp_shape
  ))
}

# Does each fitted marginal reproduce the event counts it was built from? The
# body should match almost exactly (it IS the counts); the tail is extrapolation
# and is where any disagreement will show.
for (m in list(rain_m, snow_m)) {
  if (isTRUE(m$shape_at_bound)) {
    warning(
      sprintf(
        paste0(
          "%s: GP shape is pinned at the search bound (%.3f), not an interior ",
          "optimum. Its tail ABOVE %.3f mm/h is not to be trusted; the body ",
          "below that is empirical and unaffected."
        ),
        m$col,
        m$gp_shape,
        m$u_pot
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  probe <- stats::quantile(m$events, c(0.5, 0.9, 0.97, 0.995), names = FALSE)
  message(sprintf("  %s: fitted vs empirical return period", m$col))
  for (v in probe) {
    k <- sum(m$events > v)
    message(sprintf(
      "    %7.3f mm/h -> fitted %8.2f y | empirical %8.2f y (%d events)",
      v,
      m$T_of_x(v),
      if (k > 0) m$n_years / k else NA_real_,
      k
    ))
  }
}

# %%
# ---- A grid uniform in log return period ------------------------------------
#
# Span exactly the range of raw driver values this cell's runoff peaks cover, so
# every grid point stays inside the support the forest and the copula were
# fitted on. Inverting T_of_x on a dense monotone lookup avoids ever having to
# solve it.
build_axis <- function(m, observed, n) {
  # A return period is only defined for values that qualify as declustered
  # peaks. Where a cell's runoff peaks carry driver values below the smallest
  # such peak -- near-zero rainfall at a snowmelt-driven peak, say -- those
  # points have no place on this axis and the axis starts at the smallest event
  # instead. Report how much of the cell's peak record that excludes.
  lo <- max(min(observed), m$x_min)
  hi <- max(observed)
  dropped <- mean(observed < lo)
  if (dropped > 0.005) {
    message(sprintf(
      "  note: %.1f%% of %s at this cell's runoff peaks fall below the smallest declustered event (%.3f mm/h) and are off this axis.",
      100 * dropped, m$col, m$x_min
    ))
  }
  v_dense <- seq(lo, hi, length.out = 4001)
  T_dense <- m$T_of_x(v_dense)
  keep <- c(TRUE, diff(T_dense) > 0)
  v_dense <- v_dense[keep]
  T_dense <- T_dense[keep]

  T_grid <- exp(seq(log(min(T_dense)), log(max(T_dense)), length.out = n))
  v_grid <- stats::approx(log(T_dense), v_dense, xout = log(T_grid))$y

  # dx / d(log T), for the change-of-variable Jacobian.
  dx_dlogT <- stats::approx(
    log(T_dense),
    v_dense,
    xout = log(T_grid),
    rule = 2
  )$y
  eps <- 1e-4
  up <- stats::approx(log(T_dense), v_dense, xout = log(T_grid) + eps, rule = 2)$y
  dn <- stats::approx(log(T_dense), v_dense, xout = log(T_grid) - eps, rule = 2)$y
  list(T = T_grid, x = v_grid, jac = pmax((up - dn) / (2 * eps), 0))
}

ax_rain <- build_axis(rain_m, peaks_cell$rainfall_hourly, grid_size[1])
ax_snow <- build_axis(snow_m, peaks_cell$snowmelt_hourly, grid_size[2])

message(sprintf(
  "Rain axis spans T = %.2g to %.3g y (%.2f to %.2f mm/h)",
  min(ax_rain$T), max(ax_rain$T), min(ax_rain$x), max(ax_rain$x)
))
message(sprintf(
  "Snow axis spans T = %.2g to %.3g y (%.3f to %.3f mm/h)",
  min(ax_snow$T), max(ax_snow$T), min(ax_snow$x), max(ax_snow$x)
))

gr <- expand_grid(
  ri = seq_len(grid_size[1]),
  si = seq_len(grid_size[2])
) |>
  mutate(
    rainfall_hourly = ax_rain$x[ri],
    snowmelt_hourly = ax_snow$x[si],
    T_rain = ax_rain$T[ri],
    T_snow = ax_snow$T[si],
    jacobian = ax_rain$jac[ri] * ax_snow$jac[si]
  )

# %%
frame_spec <- tibble(return_period_years = return_period_years) |>
  mutate(runoff_mm = mixture_tail_return_level(cell_mt, return_period_years * nep))

if (anyNA(frame_spec$runoff_mm)) {
  missing_T <- frame_spec$return_period_years[is.na(frame_spec$runoff_mm)]
  curve <- read_rds(
    here::here("derived", "era5_land_hourly_alps_dl_return_levels.rds")
  ) |>
    filter(cell_id == .env$cell_id, model == "GP conversion") |>
    arrange(return_period)
  frame_spec$runoff_mm[is.na(frame_spec$runoff_mm)] <- stats::approx(
    log(curve$return_period),
    curve$return_level,
    xout = log(missing_T)
  )$y
  message(sprintf(
    "%d short return period(s) filled from the cell's derived return-level curve.",
    length(missing_T)
  ))
}
stopifnot(!anyNA(frame_spec$runoff_mm))

# NB: a plain local, not a column -- pmap_dfr() below passes every column of
# frame_spec to density_at_frame() as a named argument.
extrap_ratio <- frame_spec$runoff_mm / max_obs_runoff
if (max(extrap_ratio) > 1.5) {
  warning(
    sprintf(
      paste0(
        "Cell %d: the return level reaches %.2fx the largest runoff on record ",
        "(%.2f mm/h) by T = %g y. The cell's mixture tail carries %d ",
        "component(s) with xi > 1 (max %.1f), and a mixture of GP tails is ",
        "regularly varying with index min_i(1/xi_i), so one component can set ",
        "the whole tail. Frames past that point describe the heaviest ",
        "component, not the cell. Treat them as unusable until script 4 is ",
        "rerun with gp_tail.shape_pooling: shared."
      ),
      cell_id,
      max(extrap_ratio),
      max_obs_runoff,
      max(frame_spec$return_period_years),
      sum(cell_mt$shape > 1),
      max(cell_mt$shape)
    ),
    call. = FALSE,
    immediate. = TRUE
  )
}

# %%
message("Precomputing predictive tails on the ", nrow(gr), "-point grid ...")
f_xy <- eval_joint_rain_snow_density(
  joint,
  gr$rainfall_hourly,
  gr$snowmelt_hourly
)
grid_tail <- dl_fit_cell_shared_tail(
  predict(dl_model, newdata = gr),
  shape = cell_shape
)
usable <- is.finite(grid_tail$graft_of) & is.finite(grid_tail$gp_scale)
message(sprintf("  %d of %d grid points usable", sum(usable), nrow(gr)))

conditional_density <- function(z) {
  out <- numeric(nrow(gr))
  out[usable] <- grid_tail$graft_tail_prob[usable] *
    gpd_density(
      z - grid_tail$graft_of[usable],
      grid_tail$gp_scale[usable],
      cell_shape
    )
  out
}

density_at_frame <- function(runoff_mm, return_period_years) {
  # f(rain, snow | runoff = z) times the Jacobian of (rain, snow) -> (log T_R,
  # log T_S), so the surface is a density on the plotted axes.
  num <- conditional_density(runoff_mm) * f_xy * gr$jacobian
  scale <- mean(num, na.rm = TRUE)
  tibble(
    return_period_years = return_period_years,
    runoff_mm = runoff_mm,
    extrap = runoff_mm / max_obs_runoff,
    T_rain = gr$T_rain,
    T_snow = gr$T_snow,
    rainfall_hourly = gr$rainfall_hourly,
    snowmelt_hourly = gr$snowmelt_hourly,
    density = if (is.finite(scale) && scale > 0) num / scale else NA_real_,
    density_raw = conditional_density(runoff_mm) * f_xy
  )
}

anim_tbl <- pmap_dfr(frame_spec, density_at_frame)
anim_tbl$return_period_years <- factor(
  anim_tbl$return_period_years,
  levels = frame_spec$return_period_years,
  labels = sprintf("%g", frame_spec$return_period_years)
)
# Robust shared fill breaks. pretty() over the full range collapses to a single
# band wherever one frame has a tall spike -- which is what a cell with few peak
# hours produces, its forest predictives being coarse. Resolve the bulk of the
# distribution and let one wide top band absorb the spike.
.dmax <- max(anim_tbl$density, na.rm = TRUE)
.dhi <- stats::quantile(anim_tbl$density, 0.995, na.rm = TRUE, names = FALSE)
fill_breaks <- unique(c(pretty(c(0, .dhi), n = 7), .dmax * 1.001))

# %%
# Where the mode sits, both ways, for every frame.
mode_tbl <- anim_tbl |>
  group_by(return_period_years) |>
  summarise(
    runoff_mm = first(runoff_mm),
    extrap = first(extrap),
    T_rain_mode = T_rain[which.max(density)],
    T_snow_mode = T_snow[which.max(density)],
    rain_mode = rainfall_hourly[which.max(density)],
    snow_mode = snowmelt_hourly[which.max(density)],
    T_rain_rawmode = T_rain[which.max(density_raw)],
    T_snow_rawmode = T_snow[which.max(density_raw)],
    T_rain_mean = exp(weighted.mean(log(T_rain), pmax(density, 0))),
    T_snow_mean = exp(weighted.mean(log(T_snow), pmax(density, 0))),
    .groups = "drop"
  )
print(as.data.frame(mode_tbl), digits = 3)
write_rds(
  mode_tbl,
  here::here(
    "derived",
    sprintf("rain_snow_return_period_axes_cell_%d_modes.rds", cell_id)
  )
)

# %%
rp_breaks <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500)

# Shared display window. The grid spans the full range of driver values this
# cell's runoff peaks cover -- rainfall reaches a 8700-year level -- but the
# conditional mass occupies a small part of it, so most of the panel is empty.
# Clip the VIEW (not the computation, and not per frame) to the region holding
# essentially all the mass in any frame, with a margin.
mass_window <- function(v, w, keep = 0.95, pad = 0.08) {
  o <- order(v)
  cw <- cumsum(w[o]) / sum(w[o])
  lo <- log(v[o][which.max(cw >= (1 - keep) / 2)])
  hi <- log(v[o][which.max(cw >= 1 - (1 - keep) / 2)])
  span <- hi - lo
  exp(c(lo - pad * span, hi + pad * span))
}
.w <- pmax(anim_tbl$density, 0)
.w[!is.finite(.w)] <- 0
xlim_rp <- mass_window(anim_tbl$T_rain, .w)
ylim_rp <- mass_window(anim_tbl$T_snow, .w)
message(sprintf(
  "Display window: rain T %.2g-%.3g y, snow T %.2g-%.3g y",
  xlim_rp[1], xlim_rp[2], ylim_rp[1], ylim_rp[2]
))

frame_plot <- function(rp_lab) {
  tbl <- filter(anim_tbl, return_period_years == rp_lab)
  mrow <- filter(mode_tbl, return_period_years == rp_lab)
  ggplot(tbl, aes(T_rain, T_snow)) +
    geom_contour_filled(aes(z = density), breaks = fill_breaks, alpha = 0.85) +
    geom_contour(aes(z = density), colour = "grey15", linewidth = 0.2) +
    geom_point(
      data = mrow,
      aes(T_rain_mode, T_snow_mode),
      colour = "white",
      size = 2.6,
      shape = 4,
      stroke = 1.2
    ) +
    scale_x_log10(breaks = rp_breaks, labels = \(b) ifelse(b < 1, b, as.character(b))) +
    scale_y_log10(breaks = rp_breaks, labels = \(b) ifelse(b < 1, b, as.character(b))) +
    annotation_logticks(sides = "bl", linewidth = 0.2) +
    coord_cartesian(xlim = xlim_rp, ylim = ylim_rp, expand = FALSE) +
    scale_fill_viridis_d(option = "C", end = 0.95) +
    labs(
      title = sprintf(
        "Driver return periods behind a %s-year runoff event (cell %d)",
        rp_lab,
        cell_id
      ),
      subtitle = sprintf(
        paste0(
          "runoff z = %.4g mm/h | mode: %.2g-y rainfall (%.2f mm/h) with a ",
          "%.2g-y snowmelt (%.2f mm/h)"
        ),
        mrow$runoff_mm,
        mrow$T_rain_mode,
        mrow$rain_mode,
        mrow$T_snow_mode,
        mrow$snow_mode
      ),
      x = "Rainfall return period (years, own POT)",
      y = "Snowmelt return period (years, own POT)",
      fill = NULL,
      caption = if (mrow$extrap > 1) {
        sprintf(
          paste0(
            "EXTRAPOLATION: this runoff level is %.2fx the largest ever ",
            "recorded at this cell (%.2f mm/h)."
          ),
          mrow$extrap,
          max_obs_runoff
        )
      } else {
        NULL
      }
    ) +
    theme_bw() +
    theme(plot.caption = element_text(colour = "grey25", hjust = 0))
}

message("Rendering ", nlevels(anim_tbl$return_period_years), " frames ...")
frame_paths <- map_chr(levels(anim_tbl$return_period_years), \(rp_lab) {
  path <- tempfile(pattern = "rain_snow_rp_axes_", fileext = ".png")
  ggsave(path, frame_plot(rp_lab), width = 7.2, height = 5.4, dpi = 120)
  path
})

fs::dir_create(dirname(out_gif))
frame_paths |>
  magick::image_read() |>
  magick::image_animate(fps = fps, dispose = "previous") |>
  magick::image_write(out_gif)
unlink(frame_paths)
message("Wrote ", out_gif)

# %%
message("Rendering slide stills ...")
still_paths <- map_chr(still_targets, \(target) {
  rp_lab <- sprintf("%g", target)
  stopifnot(rp_lab %in% levels(anim_tbl$return_period_years))
  path <- here::here(
    "plots",
    sprintf(
      "rain_snow_return_period_axes_cell_%d_T%03d.png",
      cell_id,
      target
    )
  )
  ggsave(
    path,
    frame_plot(rp_lab) +
      theme_bw(base_size = 15) +
      theme(plot.title = element_text(size = 16)),
    width = 12.8,
    height = 7.2,
    dpi = 200
  )
  path
})
message("Wrote:\n  ", paste(still_paths, collapse = "\n  "))
