# Joint rainfall–snowmelt explorer: marginals + copula per cell (script 6).
# Run from the package root:
#   shiny::runApp("apps/joint-rain-snow-explorer")
#
# Requires: data/era5_land_hourly_alps_all.rds and
#   data/era5_land_hourly_alps_joint_rain_snow.rds (after script 6);
#   inputs/joint_rain_snow_metadata.yaml
# Suggests: yaml, shiny, tidyverse, sf, rnaturalearth

library(shiny)
library(tidyverse)
library(sf)
library(fs)
library(rnaturalearth)
library(yaml)
library(rvinecopulib)
library(distionary)
library(famish)

repo_root <- here::here()
devtools::load_all(repo_root, quiet = TRUE)

`%||%` <- function(x, y) if (!is.null(x)) x else y

meta_path <- path(repo_root, "inputs", "joint_rain_snow_metadata.yaml")
hourly_path <- path(repo_root, "data", "era5_land_hourly_alps_all.rds")
joint_path <- path(repo_root, "data", "era5_land_hourly_alps_joint_rain_snow.rds")

margs_choices <- c(
  "empirical", "gamma", "weibull", "lnorm", "gumbel",
  "pe3", "gev", "kde"
)
bicop_family_choices <- c("parametric", "all", "onepar", "tll", "nonparametric")

default_cfg <- list(
  group_cols = list("cell_id", "x", "y"),
  rainfall_col = "rainfall_hourly",
  snowmelt_col = "snowmelt_hourly",
  marginal_rainfall = "empirical",
  marginal_snowmelt = "empirical",
  bicop_family_set = "parametric",
  bicop_controls = list(cores = 1L),
  min_obs = 40L,
  progress = FALSE,
  verbose = FALSE
)

read_joint_meta_full <- function(path) {
  if (!file.exists(path)) {
    return(list(fit_joint_rain_snow_cells = default_cfg))
  }
  yaml::read_yaml(path)
}

read_joint_cfg <- function(path) {
  meta <- read_joint_meta_full(path)
  cfg <- meta$fit_joint_rain_snow_cells
  if (is.null(cfg)) {
    return(default_cfg)
  }
  if (is.null(cfg$bicop_controls)) {
    cfg$bicop_controls <- list(cores = 1L)
  }
  cfg
}

write_joint_cfg <- function(path, cfg) {
  meta <- read_joint_meta_full(path)
  meta$fit_joint_rain_snow_cells <- cfg
  yaml::write_yaml(meta, path)
}

run_script6 <- function(root) {
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  script <- normalizePath(file.path(root, "scripts", "6-joint_rain_snow_spatial_eo.r"), mustWork = FALSE)
  if (!file.exists(script)) {
    return(list(ok = FALSE, message = paste("Script not found:", script), output = character()))
  }
  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  setwd(root)
  cmd <- paste(shQuote(normalizePath(rscript, mustWork = TRUE)), "--vanilla", shQuote(script))
  out <- tryCatch(
    suppressWarnings(system(cmd, intern = TRUE)),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    return(list(ok = FALSE, message = conditionMessage(out), output = character()))
  }
  status <- attr(out, "status")
  ok <- is.null(status) || identical(as.integer(status), 0L)
  msg_out <- paste(as.character(out), collapse = "\n")
  list(
    ok = ok,
    message = if (!ok) paste("Script exited with status", status) else "Finished.",
    output = msg_out,
    status = status
  )
}

nearest_cell <- function(lon_click, lat_click, cells_tbl) {
  dx <- cells_tbl$y - lon_click
  dy <- cells_tbl$x - lat_click
  idx <- which.min(dx^2 + dy^2)
  cells_tbl$cell_id[idx]
}

tile_dims <- function(xy_tbl, x_col = "y", y_col = "x") {
  ux <- sort(unique(xy_tbl[[x_col]]))
  uy <- sort(unique(xy_tbl[[y_col]]))
  w <- if (length(ux) > 1) stats::median(diff(ux)) else 0.25
  h <- if (length(uy) > 1) stats::median(diff(uy)) else 0.25
  c(width = w, height = h)
}

hours_per_year_from_dates <- function(dates) {
  yrs <- as.integer(format(as.Date(dates), "%Y"))
  span <- max(yrs, na.rm = TRUE) - min(yrs, na.rm = TRUE) + 1L
  span <- max(span, 1L)
  length(dates) / span
}

model_return_levels <- function(dst, hours_per_year, rp_years) {
  u <- (1 - 1 / rp_years)^(1 / hours_per_year)
  distionary::eval_quantile(dst, at = u)
}

empirical_rp_curve <- function(x, hours_per_year) {
  ok <- is.finite(x)
  x <- x[ok]
  n <- length(x)
  if (n < 2L) {
    return(tibble(x = numeric(), return_period = numeric()))
  }
  fs <- stats::ecdf(x)(x)
  fs <- pmin(pmax(fs, 1e-8), 1 - 1e-8)
  rp <- 1 / (1 - fs^hours_per_year)
  tibble(x = x, return_period = rp)
}

cdf_scores_to_copula_u <- function(u_mat, n) {
  u_mat <- pmax(pmin(as.matrix(u_mat), 1), 0)
  u_mat * (n / (n + 1L))
}

copula_density_on_normal_grid <- function(bc, ngrid = 60L) {
  if (is.null(bc)) {
    return(NULL)
  }
  eps <- 0.02
  z <- seq(qnorm(eps), qnorm(1 - eps), length.out = ngrid)
  gu <- expand.grid(z1 = z, z2 = z)
  u1 <- pnorm(gu$z1)
  u2 <- pnorm(gu$z2)
  d <- rvinecopulib::dbicop(cbind(u1, u2), bc$family, bc$rotation, bc$parameters)
  gu$d <- d
  gu
}

joint_density_original_grid <- function(joint, rain, sm, n = 45L) {
  if (is.null(joint$bicop)) {
    return(NULL)
  }
  rq <- seq(quantile(rain, 0.02, na.rm = TRUE), quantile(rain, 0.98, na.rm = TRUE), length.out = n)
  sq <- seq(quantile(sm, 0.02, na.rm = TRUE), quantile(sm, 0.98, na.rm = TRUE), length.out = n)
  gr <- expand.grid(rainfall_hourly = rq, snowmelt_hourly = sq)
  u1 <- distionary::eval_cdf(joint$marginal_rainfall, at = gr$rainfall_hourly)
  u2 <- distionary::eval_cdf(joint$marginal_snowmelt, at = gr$snowmelt_hourly)
  u1 <- pmin(pmax(as.numeric(u1), 1e-6), 1 - 1e-6)
  u2 <- pmin(pmax(as.numeric(u2), 1e-6), 1 - 1e-6)
  c_u <- rvinecopulib::dbicop(
    cbind(u1, u2),
    joint$bicop$family,
    joint$bicop$rotation,
    joint$bicop$parameters
  )
  d1 <- as.numeric(distionary::eval_density(joint$marginal_rainfall, at = gr$rainfall_hourly))
  d2 <- as.numeric(distionary::eval_density(joint$marginal_snowmelt, at = gr$snowmelt_hourly))
  gr$d <- as.numeric(c_u) * d1 * d2
  gr
}

joint0 <- read_joint_cfg(meta_path)
data_ok <- file.exists(hourly_path) && file.exists(joint_path)
hourly_all <- if (file.exists(hourly_path)) read_rds(hourly_path) else NULL
joint_tbl_init <- if (file.exists(joint_path)) read_rds(joint_path) else NULL

cells_ref <- if (!is.null(hourly_all)) {
  hourly_all |> distinct(cell_id, x, y) |> arrange(cell_id)
} else {
  tibble(cell_id = integer(), x = numeric(), y = numeric())
}

map_xlim <- if (nrow(cells_ref) > 0) range(cells_ref$y, na.rm = TRUE) + c(-0.5, 0.5) else c(-1, 1)
map_ylim <- if (nrow(cells_ref) > 0) range(cells_ref$x, na.rm = TRUE) + c(-0.5, 0.5) else c(-1, 1)
map_bbox <- st_bbox(
  c(xmin = map_xlim[1], xmax = map_xlim[2], ymin = map_ylim[1], ymax = map_ylim[2]),
  crs = st_crs(4326)
)
world_map <- tryCatch(
  rnaturalearth::ne_countries(scale = 50, returnclass = "sf") |>
    st_crop(map_bbox),
  error = function(e) NULL
)
td <- if (nrow(cells_ref) > 0) tile_dims(cells_ref) else c(width = 0.25, height = 0.25)

rp_years <- rp_reporting()

ui <- fluidPage(
  titlePanel("Joint rainfall–snowmelt explorer"),
  uiOutput("data_banner"),
  helpText(
    "Marginals are fitted with ",
    code("famish::fit_dst()"),
    "; the copula uses ",
    code("rvinecopulib::bicop()"),
    ". Edit options, write ",
    code("inputs/joint_rain_snow_metadata.yaml"),
    ", then ",
    strong("Re-run joint fitting"),
    " (",
    code("scripts/6-joint_rain_snow_spatial_eo.r"),
    "). Click the map to pick a cell."
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput(
        "marginal_rainfall",
        "Marginal: rainfall",
        choices = margs_choices,
        selected = joint0$marginal_rainfall %||% "empirical"
      ),
      selectInput(
        "marginal_snowmelt",
        "Marginal: snowmelt",
        choices = margs_choices,
        selected = joint0$marginal_snowmelt %||% "empirical"
      ),
      selectInput(
        "bicop_family_set",
        "bicop family set",
        choices = bicop_family_choices,
        selected = joint0$bicop_family_set %||% "parametric"
      ),
      numericInput("bicop_cores", "bicop cores", value = joint0$bicop_controls$cores %||% 1L, min = 1L, step = 1L),
      numericInput(
        "min_obs",
        "Minimum paired observations (copula)",
        value = joint0$min_obs %||% 40L,
        min = 10L,
        step = 1L
      ),
      radioButtons(
        "score_plane",
        "Scores / Gaussian plane",
        choices = c(
          `Rank-based normal scores (famish::nscore / implied uniform)` = "rank",
          `Marginal CDF → Gaussian (distionary::eval_cdf)` = "marginal"
        ),
        selected = "marginal"
      ),
      actionButton("rerun_joint", "Re-run joint fitting", class = "btn-primary"),
      verbatimTextOutput("rerun_log"),
      hr(),
      selectInput(
        "cell_id",
        "Cell ID",
        choices = if (nrow(cells_ref) > 0) cells_ref$cell_id else integer(),
        selected = if (nrow(cells_ref) > 0) min(cells_ref$cell_id) else NULL
      ),
      verbatimTextOutput("cell_meta")
    ),
    mainPanel(
      width = 9,
      fluidRow(
        column(6, plotOutput("map_tiles", height = "280px", click = "map_click")),
        column(6, plotOutput("fm_curve", height = "280px"))
      ),
      fluidRow(
        column(6, plotOutput("scores_plot", height = "320px")),
        column(6, plotOutput("scatter_joint", height = "320px"))
      ),
      fluidRow(
        column(6, plotOutput("hist_rain", height = "260px")),
        column(6, plotOutput("hist_snow", height = "260px"))
      )
    )
  )
)

server <- function(input, output, session) {
  joint_tbl <- reactiveVal(joint_tbl_init)
  rerun_log_txt <- reactiveVal(
    if (data_ok) {
      "Adjust settings and click Re-run to rebuild joint models."
    } else {
      "Missing hourly or joint RDS — run preprocessing scripts first."
    }
  )

  observe({
    ## Keep sidebar choices aligned when metadata exists
    cfg <- read_joint_cfg(meta_path)
    updateSelectInput(session, "marginal_rainfall", selected = cfg$marginal_rainfall %||% "empirical")
    updateSelectInput(session, "marginal_snowmelt", selected = cfg$marginal_snowmelt %||% "empirical")
    updateSelectInput(session, "bicop_family_set", selected = cfg$bicop_family_set %||% "parametric")
    updateNumericInput(session, "bicop_cores", value = cfg$bicop_controls$cores %||% 1L)
    updateNumericInput(session, "min_obs", value = cfg$min_obs %||% 40L)
  })

  output$data_banner <- renderUI({
    jt <- joint_tbl()
    if (data_ok && is.data.frame(jt) && nrow(jt) > 0L) {
      return(NULL)
    }
    missing <- c()
    if (!file.exists(hourly_path)) {
      missing <- c(missing, basename(hourly_path))
    }
    if (!file.exists(joint_path)) {
      missing <- c(missing, basename(joint_path))
    }
    div(
      class = "alert alert-warning",
      style = "padding:8px;margin-bottom:12px;",
      strong("Data not ready: "),
      paste(missing, collapse = ", "),
      " — build with tablify + joint scripts; the app still loads metadata."
    )
  })

  observeEvent(input$rerun_joint, {
    cfg <- list(
      group_cols = list("cell_id", "x", "y"),
      rainfall_col = "rainfall_hourly",
      snowmelt_col = "snowmelt_hourly",
      marginal_rainfall = input$marginal_rainfall,
      marginal_snowmelt = input$marginal_snowmelt,
      bicop_family_set = input$bicop_family_set,
      bicop_controls = list(cores = as.integer(max(1L, input$bicop_cores))),
      min_obs = as.integer(max(10L, input$min_obs)),
      progress = FALSE,
      verbose = FALSE
    )
    write_joint_cfg(meta_path, cfg)

    withProgress(message = "Running joint script…", value = 0.5, {
      res <- run_script6(repo_root)
    })

    out_txt <- res$output
    if (!nzchar(out_txt)) {
      out_txt <- res$message
    } else {
      out_txt <- paste(out_txt, res$message, sep = "\n")
    }

    if (!isTRUE(res$ok)) {
      rerun_log_txt(paste0("ERROR\n", out_txt))
      showNotification(paste("Joint script failed:", res$message), type = "error", duration = NULL)
      return(invisible(NULL))
    }

    if (file.exists(joint_path)) {
      joint_tbl(read_rds(joint_path))
    }

    rerun_log_txt(paste(
      "Last run: OK",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      paste0("(wrote ", basename(meta_path), ")"),
      "",
      out_txt,
      sep = "\n"
    ))
    showNotification("Joint fitting finished; data reloaded.", type = "message")
  })

  observeEvent(input$map_click, {
    req(nrow(cells_ref) > 0)
    cid <- nearest_cell(input$map_click$x, input$map_click$y, cells_ref)
    updateSelectInput(session, "cell_id", selected = cid)
  })

  cell_hourly <- reactive({
    req(hourly_all, input$cell_id)
    hourly_all |>
      filter(cell_id == as.integer(input$cell_id)) |>
      select(date, rainfall_hourly, snowmelt_hourly) |>
      filter(complete.cases(rainfall_hourly, snowmelt_hourly))
  })

  cell_joint <- reactive({
    req(is.data.frame(joint_tbl()))
    jt <- joint_tbl()
    ## ungroup for filter
    if (dplyr::is_grouped_df(jt)) {
      jt <- dplyr::ungroup(jt)
    }
    row <- jt |> filter(cell_id == as.integer(input$cell_id))
    if (nrow(row) != 1L) {
      return(NULL)
    }
    row$joint[[1L]]
  })

  output$cell_meta <- renderText({
    req(nrow(cells_ref) > 0)
    row <- cells_ref |> filter(cell_id == as.integer(input$cell_id))
    paste0(
      "lon = ", round(row$y, 3), ", lat = ", round(row$x, 3),
      "\n(rainfall_hourly / snowmelt_hourly in mm)"
    )
  })

  output$rerun_log <- renderText({
    rerun_log_txt()
  })

  output$map_tiles <- renderPlot({
    req(!is.null(world_map), nrow(cells_ref) > 0)
    ggplot(cells_ref, aes(y, x)) +
      geom_sf(
        data = world_map,
        inherit.aes = FALSE,
        fill = NA,
        linewidth = 1
      ) +
      geom_tile(
        aes(fill = cell_id == as.integer(input$cell_id)),
        colour = "grey30",
        linewidth = 0.3,
        width = td["width"],
        height = td["height"],
        alpha = 0.5
      ) +
      geom_text(aes(label = cell_id), size = 2.2, colour = "grey20") +
      scale_fill_manual(values = c(`TRUE` = "#b3cde3", `FALSE` = "#f7f7f7"), guide = "none") +
      coord_sf(xlim = map_xlim, ylim = map_ylim, expand = FALSE) +
      labs(x = "Longitude", y = "Latitude", title = "Cells (click to select)") +
      theme_minimal() +
      theme(panel.grid = element_blank())
  })

  output$scores_plot <- renderPlot({
    j <- cell_joint()
    h <- cell_hourly()
    validate(
      need(nrow(h) > 0, "No paired hourly data for this cell."),
      need(!is.null(j), "No joint model row for this cell — re-run script 6.")
    )
    rain <- h$rainfall_hourly
    sm <- h$snowmelt_hourly
    n <- length(rain)

    bc <- j$bicop
    validate(need(!is.null(bc), "Copula was not fitted (too few observations)."))

    if (input$score_plane == "rank") {
      z1 <- famish::nscore(rain)
      z2 <- famish::nscore(sm)
      xlab <- "Rank-based normal score (rain)"
      ylab <- "Rank-based normal score (snowmelt)"
      gu <- copula_density_on_normal_grid(bc)
    } else {
      u1 <- cdf_scores_to_copula_u(
        cbind(
          distionary::eval_cdf(j$marginal_rainfall, at = rain),
          distionary::eval_cdf(j$marginal_snowmelt, at = sm)
        ),
        n
      )
      z1 <- qnorm(u1[, 1])
      z2 <- qnorm(u1[, 2])
      xlab <- "Gaussian scores (Phi^-1(u)), rain — marginal eval_cdf"
      ylab <- "Gaussian scores (Phi^-1(u)), snow — marginal eval_cdf"
      gu <- copula_density_on_normal_grid(bc)
    }

    pt <- tibble(z1 = z1, z2 = z2)
    ggplot(pt, aes(z1, z2)) +
      geom_point(alpha = 0.35, size = 1.6, colour = "#225ea8") +
      geom_contour(
        data = gu,
        aes(z1, z2, z = d),
        inherit.aes = FALSE,
        colour = "#d94801",
        linewidth = 0.6,
        bins = 8
      ) +
      coord_cartesian(expand = FALSE) +
      theme_bw() +
      labs(
        title = "Gaussian scores vs copula density contours",
        subtitle = if (input$score_plane == "rank") {
          "Points: famish::nscore; contours: rvinecopulib::dbicop on (Phi(z1), Phi(z2))"
        } else {
          "Points and contours use marginal cdf pseudo-observations consistent with bicop()"
        },
        x = xlab,
        y = ylab
      )
  })

  output$scatter_joint <- renderPlot({
    j <- cell_joint()
    h <- cell_hourly()
    validate(
      need(nrow(h) > 0, "No paired hourly data."),
      need(!is.null(j), "No joint model row.")
    )
    rain <- h$rainfall_hourly
    sm <- h$snowmelt_hourly
    bc <- j$bicop
    validate(need(!is.null(bc), "Copula was not fitted."))

    gr <- joint_density_original_grid(j, rain, sm, n = 42L)
    ggplot(h, aes(rainfall_hourly, snowmelt_hourly)) +
      geom_point(alpha = 0.28, size = 1.2, colour = "#2c7fb8") +
      geom_contour_filled(
        data = gr,
        aes(rainfall_hourly, snowmelt_hourly, z = d),
        inherit.aes = FALSE,
        alpha = 0.55,
        bins = 10
      ) +
      scale_fill_viridis_d(option = "C", guide = guide_legend(title = "Joint\ndensity")) +
      coord_cartesian(expand = FALSE) +
      theme_bw() +
      labs(
        title = "Rainfall vs snowmelt (original scale)",
        subtitle = "Filled contours: joint density = copula density × marginal densities",
        x = "Rainfall (mm)",
        y = "Snowmelt (mm)"
      )
  })

  output$hist_rain <- renderPlot({
    j <- cell_joint()
    h <- cell_hourly()
    validate(need(nrow(h) > 0, "No data."), need(!is.null(j), "No model."))
    rng <- range(h$rainfall_hourly, na.rm = TRUE)
    xs <- seq(rng[1], max(rng[2], rng[1] + 1e-6) * 1.05, length.out = 400L)
    ys <- as.numeric(distionary::eval_density(j$marginal_rainfall, at = xs))
    hd <- tibble(x = xs, y = ys)
    ggplot(h, aes(rainfall_hourly)) +
      geom_histogram(aes(y = after_stat(density)), bins = 45, fill = "grey82", colour = "grey40") +
      geom_line(data = hd, aes(x, y), colour = "#e41a1c", linewidth = 0.9) +
      theme_bw() +
      labs(title = "Rainfall marginal", x = "Rainfall (mm)", y = "Density")
  })

  output$hist_snow <- renderPlot({
    j <- cell_joint()
    h <- cell_hourly()
    validate(need(nrow(h) > 0, "No data."), need(!is.null(j), "No model."))
    rng <- range(h$snowmelt_hourly, na.rm = TRUE)
    xs <- seq(rng[1], max(rng[2], rng[1] + 1e-6) * 1.05, length.out = 400L)
    ys <- as.numeric(distionary::eval_density(j$marginal_snowmelt, at = xs))
    hd <- tibble(x = xs, y = ys)
    ggplot(h, aes(snowmelt_hourly)) +
      geom_histogram(aes(y = after_stat(density)), bins = 45, fill = "grey82", colour = "grey40") +
      geom_line(data = hd, aes(x, y), colour = "#e41a1c", linewidth = 0.9) +
      theme_bw() +
      labs(title = "Snowmelt marginal", x = "Snowmelt (mm)", y = "Density")
  })

  output$fm_curve <- renderPlot({
    j <- cell_joint()
    h <- cell_hourly()
    validate(
      need(nrow(h) > 0, "No data."),
      need(!is.null(j), "No model.")
    )

    hpy <- hours_per_year_from_dates(h$date)

    rp <- rp_years
    mr <- model_return_levels(j$marginal_rainfall, hpy, rp)
    ms <- model_return_levels(j$marginal_snowmelt, hpy, rp)
    mod <- tibble(
      return_period = rep(rp, 2L),
      variable = rep(c("rainfall", "snowmelt"), each = length(rp)),
      level = c(mr, ms)
    )

    em_rain <- empirical_rp_curve(h$rainfall_hourly, hpy)
    em_sm <- empirical_rp_curve(h$snowmelt_hourly, hpy)
    em_rain$variable <- "rainfall"
    em_sm$variable <- "snowmelt"
    emp <- bind_rows(em_rain, em_sm) |>
      rename(level = x) |>
      filter(is.finite(return_period), return_period <= max(rp) * 3)
    if (nrow(emp) > 6000L) {
      emp <- dplyr::slice_sample(emp, n = 6000L)
    }

    ggplot(mod, aes(return_period, level, colour = variable)) +
      geom_line(linewidth = 0.9) +
      geom_point(data = emp, aes(return_period, level, colour = variable), alpha = 0.18, size = 0.9) +
      scale_x_log10(breaks = rp, minor_breaks = NULL) +
      scale_colour_manual(values = c(rainfall = "#1b9e77", snowmelt = "#d95f02")) +
      theme_bw() +
      labs(
        title = "Frequency–magnitude (marginal hourly return levels)",
        subtitle = paste0(
          "Annual recurrence (iid hours), ~",
          signif(hpy, 3),
          " h/year; empirical: ECDF-based points; model: fitted marginal quantiles"
        ),
        x = "Return period (years)",
        y = "Return level (mm)",
        colour = NULL
      ) +
      theme(legend.position = "bottom")
  })
}

shinyApp(ui, server)
