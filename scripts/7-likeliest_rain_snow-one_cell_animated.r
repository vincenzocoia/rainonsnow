# Animated contours of rainfall vs snowmelt as the runoff return period T grows.
#
# Each frame: T (years) -> z, the runoff return level from the cell's marginal,
# then a surface proportional to f(z | rain, snow) * f(rain, snow).
#
# WHY THIS IS NOW FAST
#
# The previous version evaluated, for every frame, the conditional density at
# every grid point by calling eval_density() on a distribution object: 45 x 45
# points x 30 frames is 60,000 generic calls, on top of inverting a probaverse
# mixture once per frame to get z.
#
# Both parts are closed form. Each grid point's predictive tail is a GP, so its
# density is an elementary function of (threshold, tail_prob, scale, shape);
# those are computed once, and each frame is then a single vectorised call. The
# return level z comes from the closed-form mixture tail in the same way. The
# whole animation is a few seconds, and nothing large is written to disk.
#
# The grid shares one tail shape -- the cell's, from script 4. Fitting a shape
# per grid point would reintroduce exactly the noise that
# scripts/experiments/tail-index-pooling.R measures, and would make the surface
# jump between neighbouring points for no physical reason.
#
# Requires:
#   - derived/era5_land_hourly_alps_peaks.rds (script 3)
#   - derived/era5_land_hourly_alps_dl_rqforest_models.rds (script 4)
#   - derived/era5_land_hourly_alps_dl_mixture_tails.rds (script 5)
#   - derived/era5_land_hourly_alps_joint_rain_snow.rds (script 6)
#   - magick (for GIF export)
# %%
library(tidyverse)
library(rvinecopulib)
library(yaml)
library(magick)
devtools::load_all()

# --- edit here ---
# Focus cell from inputs/rain_snow_joint_model.yaml (likeliest_rain_snow$cell_id).
# Set that to null to auto-select the most mixed rain+snowmelt cell; or override
# `cell_id` directly here.
.cfg <- read_yaml(here::here("inputs", "rain_snow_joint_model.yaml"))$likeliest_rain_snow
cell_id <- suppressWarnings(as.integer(.cfg$cell_id %||% NA_integer_))
if (is.na(cell_id)) {
  .peaks_all <- read_rds(here::here("derived", "era5_land_hourly_alps_peaks.rds"))
  .avail <- sort(unique(.peaks_all$cell_id))
  .mixed <- .peaks_all |>
    filter(rainfall_hourly > 0, snowmelt_hourly > 0) |>
    count(cell_id, sort = TRUE)
  cell_id <- if (nrow(.mixed) > 0) .mixed$cell_id[[1]] else .avail[[1]]
  message("Focus cell auto-selected (most mixed rain+snow peaks): ", cell_id)
}
# The old `runoff_cond_model` / `marginal_runoff_model` switches are gone: both
# the conditional density and the T -> z marginal now come from the cell's
# shared-shape GP tail, so there is only one coherent choice.
grid_size <- c(45L, 45L)
grid_mult <- c(1.3, 1.3)
# Slide stills are cut at these return periods. They are folded into the frame
# grid below so the stills are exact frames rather than nearest neighbours.
still_targets <- c(2, 10, 50, 200)
# Snap the nearest geometric grid point onto each target rather than appending
# it: appending would leave two frames a hair apart (exp(log(200)) is not
# exactly 200), and sprintf("%g") labels both "200", which silently merges them
# into one factor level carrying two frames' worth of points.
return_period_years <- exp(seq(log(2), log(200), length.out = 30))
for (.t in still_targets) {
  return_period_years[which.min(abs(log(return_period_years) - log(.t)))] <- .t
}
return_period_years <- sort(return_period_years)
stopifnot(!anyDuplicated(sprintf("%g", return_period_years)))
fps <- 4
out_gif <- here::here(
  "plots",
  sprintf("rain_snow_conditional_return_period_cell_%d_animated.gif", cell_id)
)
# ---

# %%
# Non-zero rain and snow only (zeros need a separate treatment later).
peaks <- read_rds(here::here("derived", "era5_land_hourly_alps_peaks.rds")) |>
  filter(rainfall_hourly != 0, snowmelt_hourly != 0)

peaks_cell <- filter(peaks, cell_id == .env$cell_id)

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
if (nrow(tail_row) != 1L) {
  stop("Expected one mixture-tail row for cell_id ", cell_id, call. = FALSE)
}
cell_mt <- tail_row$mixture_tail[[1]]
nep <- tail_row$num_events_per_year[1]

# %%
# T -> z, straight off the closed-form mixture tail.
frame_spec <- tibble(return_period_years = return_period_years) |>
  mutate(runoff_mm = mixture_tail_return_level(cell_mt, return_period_years * nep))

# Short return periods fall below the closed-form tail's base point, where the
# cell's marginal is still the empirical body of the forest mixture. Rather than
# drop those frames, read them off the return-level curve script 5 already wrote
# for this cell -- the same marginal, evaluated in its body region. Nothing is
# refitted; the curve is interpolated in log T, which is the axis it is built on.
if (anyNA(frame_spec$runoff_mm)) {
  missing_T <- frame_spec$return_period_years[is.na(frame_spec$runoff_mm)]
  curve <- read_rds(
    here::here("derived", "era5_land_hourly_alps_dl_return_levels.rds")
  ) |>
    filter(cell_id == .env$cell_id, model == "GP conversion") |>
    arrange(return_period)

  if (nrow(curve) > 1L && min(curve$return_period) <= min(missing_T)) {
    filled <- stats::approx(
      log(curve$return_period),
      curve$return_level,
      xout = log(missing_T),
      rule = 1
    )$y
    frame_spec$runoff_mm[is.na(frame_spec$runoff_mm)] <- filled
    message(sprintf(
      paste0(
        "%d short return period(s) (T <= %.1f y) sit below the closed-form ",
        "tail; filled from the cell's derived return-level curve (body region)."
      ),
      length(missing_T),
      max(missing_T)
    ))
  }

  if (anyNA(frame_spec$runoff_mm)) {
    dropped <- frame_spec$return_period_years[is.na(frame_spec$runoff_mm)]
    message(
      "Dropping ", length(dropped), " return period(s) with no return level ",
      "(T <= ", sprintf("%.1f", max(dropped)), " years)."
    )
    frame_spec <- filter(frame_spec, !is.na(runoff_mm))
  }
}
stopifnot(nrow(frame_spec) > 1L)

# %%
# Grid extent from THIS cell's peaks, not every cell's. The forest and the
# copula are both fitted on this cell alone, so a grid spanning the union of all
# cells asks them to extrapolate over a region they never saw -- and leaves most
# of the panel empty, since cell 3's peaks occupy a corner of that union.
gr <- grid_from_scatter(
  rainfall_hourly,
  snowmelt_hourly,
  data = peaks_cell,
  size = grid_size,
  mult = grid_mult
) |>
  rename(rainfall_hourly = x, snowmelt_hourly = y)

f_xy <- eval_joint_rain_snow_density(
  joint,
  gr$rainfall_hourly,
  gr$snowmelt_hourly
)

# %%
# Precompute each grid point's predictive tail ONCE, under the cell's shared
# shape. After this the per-frame work is one vectorised density evaluation.
message("Precomputing predictive tails on the ", nrow(gr), "-point grid ...")
forecast <- predict(dl_model, newdata = gr)

# ONE shape for the grid.
#
# Normally the cell's mixture already carries a single shared shape (script 4
# with gp_tail.shape_pooling: shared) and that value is simply it.
#
# The predictions file currently on disk predates that change and still carries
# a free shape per peak hour, so the mixture has one shape per component. Those
# cannot be averaged: a couple of dozen hours per cell sit above xi = 1 and a
# few above xi = 500, which is exactly the noise the shared-shape fit exists to
# remove -- it is why mixture_tail_gpd_approx(), a weighted mean, reports
# xi = 0.83 for this cell. Prefer the per-cell shared-shape estimate script 4
# wrote alongside them (shape_raw), and fall back to the weighted median of the
# components, which is robust to those outliers. The two agree closely: for
# cell 3, 0.124 against 0.130.
weighted_median <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  o <- order(x)
  x[o][which.max(cumsum(w[o]) / sum(w) >= 0.5)]
}

.shapes <- unique(round(cell_mt$shape, 10))
if (length(.shapes) == 1L) {
  cell_shape <- .shapes
  message(sprintf("Cell tail shape xi = %.3f (shared across the mixture).", cell_shape))
} else {
  .shape_file <- here::here("derived", "era5_land_hourly_alps_dl_tail_shapes.rds")
  .raw <- if (fs::file_exists(.shape_file)) {
    read_rds(.shape_file) |>
      filter(cell_id == .env$cell_id) |>
      pull(shape_raw)
  } else {
    numeric(0)
  }
  if (length(.raw) == 1L && is.finite(.raw)) {
    cell_shape <- .raw
    .src <- "script 4's per-cell shared-shape fit"
  } else {
    cell_shape <- weighted_median(cell_mt$shape, cell_mt$weights)
    .src <- "weighted median of the mixture's component shapes"
  }
  message(sprintf(
    paste0(
      "Mixture carries %d distinct component shapes (predictions were fitted ",
      "with a free shape per peak hour). Using xi = %.3f from %s."
    ),
    length(.shapes),
    cell_shape,
    .src
  ))
}
stopifnot(length(cell_shape) == 1L, is.finite(cell_shape))

grid_tail <- dl_fit_cell_shared_tail(forecast, shape = cell_shape)
usable <- is.finite(grid_tail$graft_of) & is.finite(grid_tail$gp_scale)
message(sprintf(
  "  %d of %d grid points have a usable tail (shared shape xi = %.3f)",
  sum(usable),
  nrow(gr),
  cell_shape
))

# %%
# f(z | rain, snow) for one z, over the whole grid at once.
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
  num <- conditional_density(runoff_mm) * f_xy
  # Normalising by the grid mean approximates dividing by f_Z(z), so contour
  # levels are comparable between frames.
  scale <- mean(num, na.rm = TRUE)
  tibble(
    return_period_years = return_period_years,
    runoff_mm = runoff_mm,
    rainfall_hourly = gr$rainfall_hourly,
    snowmelt_hourly = gr$snowmelt_hourly,
    density = if (is.finite(scale) && scale > 0) num / scale else NA_real_
  )
}

anim_tbl <- pmap_dfr(frame_spec, density_at_frame)
anim_tbl$return_period_years <- factor(
  anim_tbl$return_period_years,
  levels = frame_spec$return_period_years,
  labels = sprintf("%g", frame_spec$return_period_years)
)

# Shared breaks so the fill bands mean the same thing in every frame.
fill_breaks <- pretty(range(anim_tbl$density, na.rm = TRUE), n = 7)

# %%
# One plot builder, used for both the GIF frames and the standalone stills, so
# the two can never drift apart.
frame_plot <- function(rp_lab) {
  tbl <- filter(anim_tbl, return_period_years == rp_lab)
  z_mm <- unique(tbl$runoff_mm)
  ggplot(tbl, aes(rainfall_hourly, snowmelt_hourly)) +
    geom_point(data = peaks_cell, alpha = 0.15, size = 0.2) +
    geom_contour_filled(aes(z = density), breaks = fill_breaks, alpha = 0.82) +
    geom_contour(aes(z = density), colour = "grey15", linewidth = 0.2) +
    coord_cartesian(expand = FALSE) +
    scale_fill_viridis_d(option = "C", end = 0.95) +
    labs(
      title = sprintf(
        "Rainfall vs snowmelt | %s-year runoff return period (cell %d)",
        rp_lab,
        cell_id
      ),
      subtitle = sprintf(
        "z = %.4g mm/h | shared tail shape xi = %.3f",
        z_mm,
        cell_shape
      ),
      x = "Rainfall (mm/h)",
      y = "Snowmelt (mm/h)",
      fill = NULL
    ) +
    theme_bw()
}

message("Rendering ", nlevels(anim_tbl$return_period_years), " frames ...")
frame_paths <- map_chr(levels(anim_tbl$return_period_years), \(rp_lab) {
  path <- tempfile(pattern = "rain_snow_rp_", fileext = ".png")
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
# Slide stills at fixed return periods, 16:9. The animation's T grid is
# geometric, so snap each target to its nearest available frame and say so.
message("Rendering slide stills ...")
available_T <- frame_spec$return_period_years

still_paths <- map_chr(still_targets, \(target) {
  i <- which.min(abs(log(available_T) - log(target)))
  rp_lab <- levels(anim_tbl$return_period_years)[i]
  if (abs(available_T[i] / target - 1) > 0.02) {
    message(sprintf(
      "  T = %g y: nearest animation frame is T = %.3g y.",
      target,
      available_T[i]
    ))
  }
  path <- here::here(
    "plots",
    sprintf(
      "rain_snow_conditional_return_period_cell_%d_T%03d.png",
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
