# Animated contours of rainfall vs snowmelt vs runoff return period T (cell 32).
# Each frame: T (years) → z via script-5 marginal, then surface ∝ f(z|rain,snow)×f(rain,snow)
# normalized by its grid mean.
#
# Requires:
#   - derived/era5_land_hourly_alps_peaks.rds (script 3)
#   - derived/era5_land_hourly_alps_dl_rqforest_models.rds (script 4)
#   - derived/era5_land_hourly_alps_dl_marginals.rds (script 5)
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
runoff_cond_model <- "forest" # f(z | rain, snow): "forest" or "gp"
marginal_runoff_model <- "forest" # T → z marginal mixture: "forest" or "gp"
grid_size <- c(45L, 45L)
grid_mult <- c(1.3, 1.3)
# Reporting return periods T (years); z = eval_return(marginal, at = T * num_events_per_year).
return_period_years <- exp(seq(log(2), log(200), length.out = 30))
out_gif <- here::here(
  "plots",
  sprintf(
    "rain_snow_conditional_return_period_cell_%d_animated.gif",
    cell_id
  )
)
# ---

# %%
# Non-zero rain and snow only (zeros need a separate treatment later).
peaks <- read_rds(here::here("derived", "era5_land_hourly_alps_peaks.rds")) |>
  filter(rainfall_hourly != 0, snowmelt_hourly != 0)

peaks_cell <- filter(peaks, cell_id == .env$cell_id)
joint <- read_rds(here::here(
  "derived",
  "era5_land_hourly_alps_joint_rain_snow.rds"
)) |>
  filter(cell_id == .env$cell_id) |>
  pull(joint) |>
  pluck(1)

dl_model <- read_rds(
  here::here("derived", "era5_land_hourly_alps_dl_rqforest_models.rds")
) |>
  filter(cell_id == .env$cell_id) |>
  pull(dl_rqforest) |>
  pluck(1)

marg_row <- read_rds(
  here::here("derived", "era5_land_hourly_alps_dl_marginals.rds")
) |>
  filter(cell_id == .env$cell_id)
if (nrow(marg_row) != 1L) {
  stop("Expected one marginals row for cell_id ", cell_id, call. = FALSE)
}

marginal_dst <- if (marginal_runoff_model == "gp") {
  marg_row$marginal_gp[[1]]
} else {
  marg_row$marginal_forest[[1]]
}
nep <- marg_row$num_events_per_year[1]

frame_spec <- tibble(return_period_years = return_period_years) |>
  mutate(
    runoff_mm = distionary::eval_return(
      marginal_dst,
      at = return_period_years * nep
    )
  )

# %%
gr <- grid_from_scatter(
  rainfall_hourly,
  snowmelt_hourly,
  data = peaks,
  size = grid_size,
  mult = grid_mult
) |>
  rename(rainfall_hourly = x, snowmelt_hourly = y)

f_xy <- eval_joint_rain_snow_density(
  joint,
  gr$rainfall_hourly,
  gr$snowmelt_hourly
)

forecast <- predict(dl_model, newdata = gr)
if (runoff_cond_model == "forest") {
  forecast <- purrr::map(
    forecast,
    \(d) tryCatch(convert_emp_to_gp(d), error = function(e) d),
    .progress = TRUE
  )
}

# %%
# Per-z numerator and mean normalization (≈ f(rain, snow | z) up to a z-dependent
# constant). On a fine grid, mean(f(z|xy)*f(xy)) tracks ∫∫ f(z|xy)f(xy) dx dy,
# i.e. f_Z(z) under the mixture/grid approximation—so dividing by the mean is a
# practical substitute for loading marginal f_Z from script 5.

density_at_frame <- function(runoff_mm, return_period_years) {
  f_z_given_xy <- purrr::map_dbl(forecast, \(d) {
    v <- distionary::eval_density(d, at = runoff_mm)
    if (length(v) != 1L || !is.finite(v)) {
      return(0)
    }
    v
  })
  num <- f_z_given_xy * f_xy
  tibble(
    return_period_years = return_period_years,
    runoff_mm = runoff_mm,
    rainfall_hourly = gr$rainfall_hourly,
    snowmelt_hourly = gr$snowmelt_hourly,
    density = num / mean(num)
  )
}

anim_tbl <- purrr::pmap_dfr(frame_spec, density_at_frame, .progress = TRUE)
anim_tbl$return_period_years <- factor(
  anim_tbl$return_period_years,
  levels = frame_spec$return_period_years,
  labels = sprintf("%g", frame_spec$return_period_years)
)

# Shared breaks so fill bands match across frames (geom_contour_filled is discrete).
fill_breaks <- pretty(range(anim_tbl$density, na.rm = TRUE), n = 7)

# %%
frame_paths <- purrr::map_chr(levels(anim_tbl$return_period_years), \(rp_lab) {
  tbl <- filter(anim_tbl, return_period_years == rp_lab)
  z_mm <- unique(tbl$runoff_mm)
  p <- ggplot(tbl, aes(rainfall_hourly, snowmelt_hourly)) +
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
        "z = %.4g mm/h | f(rain, snow | z) ≈ f(z|rain,snow)×f(rain,snow) / grid mean",
        z_mm
      ),
      x = "Rainfall (mm/h)",
      y = "Snowmelt (mm/h)",
      fill = NULL
    ) +
    theme_bw()

  path <- tempfile(pattern = "rain_snow_rp_", fileext = ".png")
  ggsave(path, p, width = 7.2, height = 5.4, dpi = 120)
  path
})

fs::dir_create(dirname(out_gif), showWarnings = FALSE, recursive = TRUE)
frame_paths |>
  magick::image_read() |>
  magick::image_animate(fps = 4, dispose = "previous") |>
  magick::image_write(out_gif)

unlink(frame_paths)
message("Wrote ", out_gif)
