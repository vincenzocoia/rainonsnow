# Conditional density of hourly rainfall vs snowmelt given hourly runoff = z:
# proportional to f(z | rain, snow) × f(rain, snow), using the rqforest from
# script 4 and the joint marginal+copula model from script 6.
#
# Requires:
#   - inputs/rain_snow_joint_model.yaml (likeliest_rain_snow block)
#   - derived/era5_land_hourly_alps_all.rds (script 2)
#   - derived/era5_land_hourly_alps_dl_rqforest_models.rds (script 4)
#   - derived/era5_land_hourly_alps_joint_rain_snow.rds (script 6)
# Optional for return-period choice of z:
#   - derived/era5_land_hourly_alps_dl_return_levels.rds (script 5)
# Optional for normalization by f_Z(z):
#   - derived/era5_land_hourly_alps_dl_marginals.rds (script 5)
# %%
library(tidyverse)
library(rlang)
library(yaml)
library(rvinecopulib)
library(logger)
devtools::load_all()

# %%
meta_path <- here::here("inputs", "rain_snow_joint_model.yaml")
if (!file.exists(meta_path)) {
  stop("Missing ", meta_path, "; see repository template.", call. = FALSE)
}

meta <- read_yaml(meta_path)
cfg <- meta$likeliest_rain_snow
if (is.null(cfg) || !is.list(cfg)) {
  stop(
    "Add a `likeliest_rain_snow:` block to inputs/rain_snow_joint_model.yaml ",
    "(see comments in that file).",
    call. = FALSE
  )
}

cell_id <- as.integer(cfg$cell_id %||% NA_integer_)
if (is.na(cell_id)) {
  stop(
    "`cell_id` must be set under likeliest_rain_snow in inputs/rain_snow_joint_model.yaml",
    call. = FALSE
  )
}

marginal_rp_model <- match.arg(
  tolower(as.character(cfg$marginal_return_model %||% "gp")),
  c("gp", "forest")
)
runoff_cond_model <- match.arg(
  tolower(as.character(cfg$runoff_conditional_model %||% "forest")),
  c("forest", "gp")
)
normalize_fz <- isTRUE(cfg$normalize_by_marginal_runoff_density)

grid_size <- as.integer(cfg$grid$size %||% c(45L, 45L))
grid_mult <- cfg$grid$mult %||% c(1.3, 1.3)

runoff_mm <- cfg$runoff_threshold_mm
rp_years <- cfg$return_period

# %%
# Target runoff z (mm/h): fixed threshold or return level from script 5.

if (
  !is.null(runoff_mm) &&
    length(runoff_mm) == 1L &&
    is.finite(as.numeric(runoff_mm))
) {
  z <- as.numeric(runoff_mm)
  z_src <- sprintf("fixed threshold (%g mm/h)", z)
} else if (!is.null(rp_years) && is.finite(as.numeric(rp_years))) {
  rp_years <- as.numeric(rp_years)
  lv_path <- here::here(
    "derived",
    "era5_land_hourly_alps_dl_return_levels.rds"
  )
  if (!file.exists(lv_path)) {
    stop(
      "Need script 5 output ",
      basename(lv_path),
      " to resolve return_period, or set `runoff_threshold_mm` instead.",
      call. = FALSE
    )
  }

  marginal_rp_label <- switch(
    marginal_rp_model,
    gp = "GP conversion",
    forest = "Random Forest"
  )

  lvls <- read_rds(lv_path)
  if (
    "return_period" %in% names(lvls) && !"return_period_years" %in% names(lvls)
  ) {
    lvls <- dplyr::rename(lvls, return_period_years = return_period)
  }

  mag <- lvls |>
    dplyr::filter(
      .data$cell_id == .env$cell_id,
      abs(.data$return_period_years - .env$rp_years) < 1e-4,
      .data$model == .env$marginal_rp_label
    )
  if (nrow(mag) != 1L) {
    stop(
      "Could not find return level for cell ",
      cell_id,
      ", return_period ",
      rp_years,
      ", model ",
      marginal_rp_label,
      call. = FALSE
    )
  }

  z <- mag$return_level[1]
  z_src <- sprintf(
    "%g-year return level (%s; %g mm/h)",
    rp_years,
    marginal_rp_label,
    z
  )
} else {
  stop(
    "Set either `runoff_threshold_mm` or `return_period` under likeliest_rain_snow in inputs/rain_snow_joint_model.yaml.",
    call. = FALSE
  )
}

log_info("Starting 7-likeliest_rain_snow.r")
log_info(paste("Cell", cell_id, "| z =", signif(z, 6), "mm/h —", z_src))

# %%
# Peaks table, joint rain–snow model, and rqforest for f(z | rain, snow).

peaks_path <- here::here("derived", "era5_land_hourly_alps_peaks.rds")
joint_path <- here::here("derived", "era5_land_hourly_alps_joint_rain_snow.rds")
models_path <- here::here(
  "derived",
  "era5_land_hourly_alps_dl_rqforest_models.rds"
)

for (p in c(peaks_path, joint_path, models_path)) {
  if (!file.exists(p)) {
    stop("Missing required file: ", p, call. = FALSE)
  }
}

peaks <- read_rds(peaks_path)
joint_tbl <- read_rds(joint_path)
models_tbl <- read_rds(models_path)

peaks <- filter(peaks, rainfall_hourly != 0, snowmelt_hourly != 0)

peaks_cell <- dplyr::filter(peaks, .data$cell_id == .env$cell_id)
if (nrow(peaks_cell) == 0L) {
  stop("No hourly rows for cell_id ", cell_id, call. = FALSE)
}

joint_row <- dplyr::filter(joint_tbl, .data$cell_id == .env$cell_id)
if (nrow(joint_row) != 1L) {
  stop(
    "Expected exactly one joint model row for cell_id ",
    cell_id,
    call. = FALSE
  )
}
joint <- joint_row$joint[[1]]

model_row <- dplyr::filter(models_tbl, .data$cell_id == .env$cell_id)
if (nrow(model_row) != 1L) {
  stop(
    "Expected exactly one dl_rqforest row for cell_id ",
    cell_id,
    call. = FALSE
  )
}
dl_model <- model_row$dl_rqforest[[1]]

# %%
# Regular grid over the observed rain–snow range for this cell.

gr <- grid_from_scatter(
  rainfall_hourly,
  snowmelt_hourly,
  data = peaks,
  size = grid_size,
  mult = grid_mult
) |>
  rename(rainfall_hourly = x, snowmelt_hourly = y)

# %%
# Joint density f(rain, snow) from script 6.

f_xy <- eval_joint_rain_snow_density(
  joint,
  gr$rainfall_hourly,
  gr$snowmelt_hourly
)

gr |>
  mutate(joint_density = f_xy) |>
  ggplot(aes(rainfall_hourly, snowmelt_hourly)) +
  geom_point(data = peaks, alpha = 0.2, size = 0.2) +
  geom_contour_filled(aes(z = log(joint_density)), alpha = 0.82) +
  theme_bw()
# %%
# Conditional runoff density f(z | rain, snow) at each grid point (script 4).

forecast <- predict(dl_model, newdata = gr)
if (runoff_cond_model == "forest") {
  forecast <- purrr::map(
    forecast,
    \(d) tryCatch(convert_emp_to_gp(d), error = function(e) d),
    .progress = TRUE
  )
}

f_z_given_xy <- purrr::map_dbl(forecast, \(d) {
  v <- distionary::eval_density(d, at = z)
  if (length(v) != 1L || !is.finite(v)) {
    return(0)
  }
  v
})

# Take a look at f_z_given_xy:
gr |>
  mutate(f_z_given_xy = f_z_given_xy) |>
  ggplot(aes(rainfall_hourly, snowmelt_hourly)) +
  geom_point(data = peaks, alpha = 0.2, size = 0.2) +
  geom_contour_filled(aes(z = (f_z_given_xy)), alpha = 0.82) +
  theme_bw()

# %%
# Surface proportional to f(rain, snow | z); optionally divide by f_Z(z).

surface <- f_z_given_xy * f_xy
marginal_density_z <- NA_real_

if (normalize_fz) {
  marg_path <- here::here("derived", "era5_land_hourly_alps_dl_marginals.rds")
  if (!file.exists(marg_path)) {
    stop(
      "`normalize_by_marginal_runoff_density: true` requires ",
      basename(marg_path),
      " from script 5.",
      call. = FALSE
    )
  }

  marg_row <- read_rds(marg_path) |>
    dplyr::filter(.data$cell_id == .env$cell_id)
  if (nrow(marg_row) != 1L) {
    stop("Marginal mixtures: no single row for cell ", cell_id, call. = FALSE)
  }

  marginal_dst <- if (runoff_cond_model == "gp") {
    marg_row$marginal_gp[[1]]
  } else {
    marg_row$marginal_forest[[1]]
  }

  marginal_density_z <- as.numeric(distionary::eval_density(
    marginal_dst,
    at = z
  ))

  if (!is.finite(marginal_density_z) || marginal_density_z <= 0) {
    warning(
      "Marginal density f_Z(z) is not usable; plotting unnormalized surface.",
      call. = FALSE
    )
  } else {
    surface <- surface / marginal_density_z
  }
}

plot_tbl <- mutate(
  gr,
  joint_xy = f_xy,
  fz_given_xy = f_z_given_xy,
  density_propto = surface
)

# %%
density_label <- if (
  normalize_fz && is.finite(marginal_density_z) && marginal_density_z > 0
) {
  sprintf("f(rain, snow | runoff = %.4g)", z)
} else {
  sprintf(
    "density proportional to f(z|rain,snow) * f(rain,snow), z = %.4g mm/h",
    z
  )
}

p <- ggplot(plot_tbl, aes(rainfall_hourly, snowmelt_hourly)) +
  geom_contour_filled(aes(z = density_propto), alpha = 0.82) +
  geom_contour(
    aes(z = density_propto),
    colour = "grey15",
    linewidth = 0.2
  ) +
  coord_cartesian(expand = FALSE) +
  scale_fill_viridis_d(option = "C", end = 0.95) +
  labs(
    title = paste0(
      "Rainfall vs snowmelt | hourly runoff @ ",
      signif(z, 4),
      " mm/h"
    ),
    subtitle = paste0(density_label, "\n", z_src),
    x = "Rainfall (mm/h)",
    y = "Snowmelt (mm/h)",
    fill = NULL
  ) +
  theme_bw()

# %%
plots_dir <- here::here("plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

slug_z <- gsub("[^0-9eE+.-]", "_", sprintf("%.4g", z))

out_pdf <- cfg$output_pdf %||%
  file.path(
    plots_dir,
    sprintf("rain_snow_conditional_runoff_cell_%d_z_%s.pdf", cell_id, slug_z)
  )
ggplot2::ggsave(out_pdf, p, width = 7.2, height = 5.4)
log_info(paste("Wrote", out_pdf))

out_rds <- cfg$output_rds %||%
  here::here(
    "derived",
    sprintf("rain_snow_conditional_runoff_cell_%d_z_%s.rds", cell_id, slug_z)
  )

saveRDS(
  list(
    cell_id = cell_id,
    runoff_mm = z,
    z_source = z_src,
    joint_model_note = if (is.null(joint$bicop)) {
      "independence (no copula fitted)"
    } else {
      "bicop joint density"
    },
    runoff_conditional_model = runoff_cond_model,
    normalized = normalize_fz,
    marginal_fz = marginal_density_z,
    grid = plot_tbl
  ),
  out_rds
)
log_info(paste("Wrote", out_rds))
log_info("Finished 7-likeliest_rain_snow.r")
