# Conditional density of hourly rainfall vs snowmelt given hourly runoff = z:
# proportional to f(z | rain, snow) × f(rain, snow), using the rqforest from
# script 4 and the joint marginal+copula model from script 6.
#
# Requires:
#   - data/era5_land_hourly_alps_all.rds (script 2)
#   - data/era5_land_hourly_alps_dl_rqforest_models.rds (script 4)
#   - data/era5_land_hourly_alps_joint_rain_snow.rds (script 6)
# Optional for return-period choice of z:
#   - data/era5_land_hourly_alps_dl_marginal_return_levels.rds (script 5)
# Optional for normalization by f_Z(z):
#   - data/era5_land_hourly_alps_dl_marginals.rds (script 5)
# %%
library(tidyverse)
library(rlang)
library(yaml)
library(rvinecopulib)
library(logger)
devtools::load_all()

meta_path <- here::here("inputs", "rain_snow_conditional_runoff.yaml")
if (!file.exists(meta_path)) {
  stop("Missing ", meta_path, "; see repository template.", call. = FALSE)
}

meta <- read_yaml(meta_path)
cfg <- meta$rain_snow_conditional_runoff

cell_id <- as.integer(cfg$cell_id %||% NA_integer_)
if (is.na(cell_id)) {
  stop("`cell_id` must be set in inputs/rain_snow_conditional_runoff.yaml", call. = FALSE)
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

grid_size <- cfg$grid$size %||% c(45L, 45L)
grid_mult <- cfg$grid$mult %||% c(1.3, 1.3)
grid_size <- as.integer(grid_size)

runoff_mm <- cfg$runoff_threshold_mm
rp_years <- cfg$return_period

if (!is.null(runoff_mm) && length(runoff_mm) == 1L && is.finite(as.numeric(runoff_mm))) {
  z <- as.numeric(runoff_mm)
  z_src <- sprintf("fixed threshold (%g mm/h)", z)
} else if (!is.null(rp_years) && is.finite(as.numeric(rp_years))) {
  rp_years <- as.numeric(rp_years)
  lv_path <- here::here("data", "era5_land_hourly_alps_dl_marginal_return_levels.rds")
  if (!file.exists(lv_path)) {
    stop(
      "Need script 5 output ",
      basename(lv_path),
      " to resolve return_period, or set `runoff_threshold_mm` instead.",
      call. = FALSE
    )
  }
  bundle <- read_rds(lv_path)
  lvls <- bundle$marginal_return_levels_long
  model_lab <- switch(
    marginal_rp_model,
    gp = "Forest mixture + GP tail",
    forest = "Forest mixture (empirical)"
  )
  mag <- lvls |>
    dplyr::filter(
      .data$cell_id == .env$cell_id,
      abs(.data$return_period_years - .env$rp_years) < 1e-4,
      .data$model == .env$model_lab
    )
  if (nrow(mag) != 1L) {
    stop(
      "Could not find return level for cell ", cell_id,
      ", return_period ", rp_years, ", model ", model_lab,
      call. = FALSE
    )
  }
  z <- mag$return_level[1]
  z_src <- sprintf("%g-year return level (%s; %g mm/h)", rp_years, model_lab, z)
} else {
  stop(
    "Set either `runoff_threshold_mm` or `return_period` in the yaml.",
    call. = FALSE
  )
}

log_info("Starting 7-likeliest_rain_snow.r")
log_info(paste("Cell", cell_id, "| z =", signif(z, 6), "mm/h —", z_src))

hourly_path <- here::here("data", "era5_land_hourly_alps_all.rds")
joint_path <- here::here("data", "era5_land_hourly_alps_joint_rain_snow.rds")
models_path <- here::here("data", "era5_land_hourly_alps_dl_rqforest_models.rds")

for (p in c(hourly_path, joint_path, models_path)) {
  if (!file.exists(p)) {
    stop("Missing required file: ", p, call. = FALSE)
  }
}

hourly_all <- read_rds(hourly_path)
joint_tbl <- read_rds(joint_path)
models_tbl <- read_rds(models_path)

hourly_cell <- hourly_all |> dplyr::filter(.data$cell_id == .env$cell_id)
if (nrow(hourly_cell) == 0L) {
  stop("No hourly rows for cell_id ", cell_id, call. = FALSE)
}

joint_row <- joint_tbl |> dplyr::filter(.data$cell_id == .env$cell_id)
if (nrow(joint_row) != 1L) {
  stop("Expected exactly one joint model row for cell_id ", cell_id, call. = FALSE)
}
joint <- joint_row$joint[[1]]

model_row <- models_tbl |> dplyr::filter(.data$cell_id == .env$cell_id)
if (nrow(model_row) != 1L) {
  stop("Expected exactly one dl_rqforest row for cell_id ", cell_id, call. = FALSE)
}
dl_model <- model_row$dl_rqforest[[1]]

gr <- grid_from_scatter(
  rainfall_hourly,
  snowmelt_hourly,
  data = hourly_cell,
  size = grid_size,
  mult = grid_mult
) |>
  dplyr::rename(rainfall_hourly = x, snowmelt_hourly = y)

eval_joint_density_xy <- function(joint, rainfall_hourly, snowmelt_hourly) {
  if (is.null(joint$bicop)) {
    d1 <- as.numeric(distionary::eval_density(joint$marginal_rainfall, at = rainfall_hourly))
    d2 <- as.numeric(distionary::eval_density(joint$marginal_snowmelt, at = snowmelt_hourly))
    return(d1 * d2)
  }
  u1 <- distionary::eval_cdf(joint$marginal_rainfall, at = rainfall_hourly)
  u2 <- distionary::eval_cdf(joint$marginal_snowmelt, at = snowmelt_hourly)
  u1 <- pmin(pmax(as.numeric(u1), 1e-6), 1 - 1e-6)
  u2 <- pmin(pmax(as.numeric(u2), 1e-6), 1 - 1e-6)
  c_u <- rvinecopulib::dbicop(
    cbind(u1, u2),
    joint$bicop$family,
    joint$bicop$rotation,
    joint$bicop$parameters
  )
  d1 <- as.numeric(distionary::eval_density(joint$marginal_rainfall, at = rainfall_hourly))
  d2 <- as.numeric(distionary::eval_density(joint$marginal_snowmelt, at = snowmelt_hourly))
  as.numeric(c_u) * d1 * d2
}

f_xy <- eval_joint_density_xy(joint, gr$rainfall_hourly, gr$snowmelt_hourly)

forecast <- predict(dl_model, newdata = gr)
if (runoff_cond_model == "gp") {
  forecast <- purrr::map(
    forecast,
    function(d) {
      tryCatch(
        convert_emp_to_gp(d),
        error = function(e) d
      )
    }
  )
}

f_z_given_xy <- purrr::map_dbl(forecast, function(d) {
  v <- suppressWarnings(as.numeric(distionary::eval_density(d, at = z)))
  if (length(v) != 1L || !is.finite(v)) {
    return(0)
  }
  v
})

surface <- f_z_given_xy * f_xy

marginal_density_z <- NA_real_
if (normalize_fz) {
  marg_path <- here::here("data", "era5_land_hourly_alps_dl_marginals.rds")
  if (!file.exists(marg_path)) {
    stop(
      "`normalize_by_marginal_runoff_density: true` requires ",
      basename(marg_path),
      " from script 5.",
      call. = FALSE
    )
  }
  marg_tbl <- read_rds(marg_path)
  mrow <- marg_tbl |> dplyr::filter(.data$cell_id == .env$cell_id)
  if (nrow(mrow) != 1L) {
    stop("Marginal mixtures: no single row for cell ", cell_id, call. = FALSE)
  }
  marginal_dst <- if (runoff_cond_model == "gp") {
    mrow$marginal_gp[[1]]
  } else {
    mrow$marginal_forest[[1]]
  }
  marginal_density_z <- as.numeric(distionary::eval_density(marginal_dst, at = z))
  if (!is.finite(marginal_density_z) || marginal_density_z <= 0) {
    warning(
      "Marginal density f_Z(z) is not usable; plotting unnormalized surface.",
      call. = FALSE
    )
  } else {
    surface <- surface / marginal_density_z
  }
}

plot_tbl <- gr |>
  dplyr::mutate(
    joint_xy = f_xy,
    fz_given_xy = f_z_given_xy,
    density_propto = surface
  )

lab_z <- if (normalize_fz && is.finite(marginal_density_z) && marginal_density_z > 0) {
  sprintf("f(rain, snow | runoff = %.4g)", z)
} else {
  sprintf(
    "density proportional to f(z|rain,snow) * f(rain,snow), z = %.4g mm/h",
    z
  )
}

p <- ggplot(plot_tbl, aes(.data$rainfall_hourly, .data$snowmelt_hourly)) +
  geom_contour_filled(aes(z = .data$density_propto), alpha = 0.82) +
  geom_contour(aes(z = .data$density_propto), colour = alpha("grey15", 0.35), linewidth = 0.2) +
  coord_cartesian(expand = FALSE) +
  scale_fill_viridis_d(option = "C", end = 0.95) +
  labs(
    title = paste0("Rainfall vs snowmelt | hourly runoff @ ", signif(z, 4), " mm/h"),
    subtitle = paste0(lab_z, "\n", z_src),
    x = "Rainfall (mm/h)",
    y = "Snowmelt (mm/h)",
    fill = NULL
  ) +
  theme_bw()

plots_dir <- here::here("plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

slug_z <- sprintf("%.4g", z)
slug_z <- gsub("[^0-9eE+.-]", "_", slug_z)

out_pdf <- cfg$output_pdf %||% file.path(
  plots_dir,
  sprintf("rain_snow_conditional_runoff_cell_%d_z_%s.pdf", cell_id, slug_z)
)
ggplot2::ggsave(out_pdf, p, width = 7.2, height = 5.4)
log_info(paste("Wrote", out_pdf))

out_rds <- cfg$output_rds %||% here::here(
  "data",
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
