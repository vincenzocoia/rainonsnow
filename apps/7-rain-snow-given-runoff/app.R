# Rain-snow structure given an extreme runoff event.
#
# The interactive version of scripts/7-likeliest_rain_snow-one_cell_animated.r.
# Instead of rendering a GIF, the return period is a slider with a play button,
# so you can scrub T as well as watch it run.
#
# Run from the package root:
#   shiny::runApp("apps/7-rain-snow-given-runoff")
#
# Requires:
#   - derived/era5_land_hourly_alps_peaks.rds            (script 3)
#   - derived/era5_land_hourly_alps_dl_rqforest_models.rds (script 4)
#   - derived/era5_land_hourly_alps_dl_mixture_tails.rds  (script 5)
#   - derived/era5_land_hourly_alps_joint_rain_snow.rds   (script 6)
#
# This is only interactive because the tail machinery is closed form: each grid
# point's predictive tail is fitted ONCE per cell, and every frame is then a
# single vectorised density evaluation. See R/mixture_tail.R.
#
# Suggests: shiny, tidyverse, fs, rvinecopulib

library(shiny)
library(tidyverse)
library(fs)
library(rvinecopulib)

repo_root <- here::here()
devtools::load_all(repo_root, quiet = TRUE)
source(file.path(repo_root, "apps", "ros_theme.R"))

peaks_path <- path(repo_root, "derived", "era5_land_hourly_alps_peaks.rds")
models_path <- path(repo_root, "derived", "era5_land_hourly_alps_dl_rqforest_models.rds")
tails_path <- path(repo_root, "derived", "era5_land_hourly_alps_dl_mixture_tails.rds")
joint_path <- path(repo_root, "derived", "era5_land_hourly_alps_joint_rain_snow.rds")

read_or_null <- function(p) if (file.exists(p)) readRDS(p) else NULL

peaks_all <- read_or_null(peaks_path)
models_tbl <- read_or_null(models_path)
tails_tbl <- read_or_null(tails_path)
joint_tbl <- read_or_null(joint_path)

missing_files <- c(
  if (is.null(peaks_all)) basename(peaks_path),
  if (is.null(models_tbl)) basename(models_path),
  if (is.null(tails_tbl)) basename(tails_path),
  if (is.null(joint_tbl)) basename(joint_path)
)

peaks_nz <- if (!is.null(peaks_all)) {
  dplyr::filter(peaks_all, rainfall_hourly != 0, snowmelt_hourly != 0)
} else {
  NULL
}

cell_choices <- if (!is.null(tails_tbl)) sort(unique(tails_tbl$cell_id)) else integer()

# Default to the cell with the most co-occurring rain and snowmelt peaks, the
# same rule scripts 6 and 7 use when likeliest_rain_snow$cell_id is null.
default_cell <- if (!is.null(peaks_nz) && nrow(peaks_nz) > 0) {
  counts <- dplyr::count(peaks_nz, cell_id, sort = TRUE)
  candidates <- intersect(counts$cell_id, cell_choices)
  if (length(candidates)) candidates[[1]] else NULL
} else if (length(cell_choices)) {
  cell_choices[[1]]
} else {
  NULL
}

ui <- fluidPage(
  theme = ros_bs_theme(),
  ros_header(
    "Rain and snowmelt given extreme runoff",
    "How the likely combination of rainfall and snowmelt shifts as the runoff return period grows."
  ),
  uiOutput("banner"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput(
        "cell_id",
        "Cell ID",
        choices = cell_choices,
        selected = default_cell
      ),
      sliderInput(
        "log_rp",
        "Return period T (years)",
        min = log10(2),
        max = log10(200),
        value = log10(10),
        step = (log10(200) - log10(2)) / 29,
        ticks = FALSE,
        animate = animationOptions(interval = 450, loop = TRUE)
      ),
      helpText(
        "Press play to sweep T. The slider is on a log scale; the value in the ",
        "subtitle is the runoff return level it corresponds to."
      ),
      sliderInput("grid_n", "Grid resolution", min = 25, max = 75, value = 45, step = 5),
      checkboxInput("show_points", "Show this cell's peak hours", TRUE),
      ros_note(
        "The surface is proportional to f(runoff | rain, snow) x f(rain, snow), ",
        "normalised by its own grid mean so the contour bands mean the same ",
        "thing at every T. The conditional density uses the cell's shared tail ",
        "shape, so the surface does not jump between neighbouring grid points."
      ),
      verbatimTextOutput("cell_meta")
    ),
    mainPanel(
      width = 9,
      plotOutput("surface", height = "560px"),
      plotOutput("return_curve", height = "230px")
    )
  )
)

server <- function(input, output, session) {
  output$banner <- renderUI({
    if (!length(missing_files)) {
      return(NULL)
    }
    div(
      class = "ros-note",
      "Missing derived files: ",
      paste(missing_files, collapse = ", "),
      ". Run the pipeline through script 6."
    )
  })

  cell_bits <- reactive({
    req(length(missing_files) == 0, input$cell_id)
    cid <- as.integer(input$cell_id)

    trow <- dplyr::filter(tails_tbl, cell_id == cid)
    mrow <- dplyr::filter(models_tbl, cell_id == cid)
    jrow <- dplyr::filter(joint_tbl, cell_id == cid)
    validate(
      need(nrow(trow) == 1L, "No mixture tail for this cell — rerun script 5."),
      need(nrow(mrow) == 1L, "No fitted model for this cell — rerun script 4."),
      need(nrow(jrow) == 1L, "No joint model for this cell — rerun script 6.")
    )
    list(
      mt = trow$mixture_tail[[1]],
      nep = trow$num_events_per_year[1],
      model = mrow$dl_rqforest[[1]],
      joint = jrow$joint[[1]],
      peaks = dplyr::filter(peaks_nz, cell_id == cid)
    )
  })

  # Everything that does not depend on T is computed once per cell (and per grid
  # size), so moving the slider only re-evaluates a vectorised density.
  surface_base <- reactive({
    b <- cell_bits()
    n <- as.integer(input$grid_n)

    gr <- grid_from_scatter(
      rainfall_hourly,
      snowmelt_hourly,
      data = peaks_nz,
      size = c(n, n),
      mult = c(1.3, 1.3)
    ) |>
      dplyr::rename(rainfall_hourly = x, snowmelt_hourly = y)

    f_xy <- eval_joint_rain_snow_density(
      b$joint,
      gr$rainfall_hourly,
      gr$snowmelt_hourly
    )

    shape <- unique(round(b$mt$shape, 10))
    validate(need(
      length(shape) == 1L,
      "This cell's tail shapes are not shared; set gp_tail.shape_pooling: shared and rerun script 4."
    ))

    withProgress(message = "Fitting predictive tails on the grid", value = 0.5, {
      forecast <- predict(b$model, newdata = gr)
      gt <- dl_fit_cell_shared_tail(forecast, shape = shape)
    })
    usable <- is.finite(gt$graft_of) & is.finite(gt$gp_scale)

    list(gr = gr, f_xy = f_xy, gt = gt, usable = usable, shape = shape, bits = b)
  })

  runoff_level <- reactive({
    b <- cell_bits()
    rp <- 10^input$log_rp
    z <- mixture_tail_return_level(b$mt, rp * b$nep)
    list(rp = rp, z = z)
  })

  output$surface <- renderPlot({
    s <- surface_base()
    rl <- runoff_level()
    validate(need(
      is.finite(rl$z),
      sprintf(
        "T = %.1f years is below this cell's tail region (valid for T above about %.1f years).",
        rl$rp,
        1 / (mixture_tail_max_prob(s$bits$mt) * s$bits$nep)
      )
    ))

    dens <- numeric(nrow(s$gr))
    dens[s$usable] <- s$gt$graft_tail_prob[s$usable] *
      gpd_density(
        rl$z - s$gt$graft_of[s$usable],
        s$gt$gp_scale[s$usable],
        s$shape
      )
    num <- dens * s$f_xy
    scale_by <- mean(num, na.rm = TRUE)
    validate(need(
      is.finite(scale_by) && scale_by > 0,
      "The surface is numerically zero everywhere at this return period."
    ))

    tbl <- s$gr
    tbl$density <- num / scale_by

    p <- ggplot(tbl, aes(rainfall_hourly, snowmelt_hourly))
    if (isTRUE(input$show_points)) {
      p <- p + geom_point(
        data = s$bits$peaks,
        alpha = 0.16,
        size = 0.5,
        colour = ros_palette$ink
      )
    }
    p +
      geom_contour_filled(aes(z = density), bins = 9, alpha = 0.8) +
      geom_contour(aes(z = density), colour = "grey20", linewidth = 0.2, bins = 9) +
      scale_fill_viridis_d(option = "C", end = 0.95) +
      coord_cartesian(expand = FALSE) +
      labs(
        title = sprintf("T = %.1f-year runoff, cell %s", rl$rp, input$cell_id),
        subtitle = sprintf(
          "runoff level %.4g mm/h  |  shared tail shape xi = %.3f",
          rl$z,
          s$shape
        ),
        x = "Rainfall (mm/h)",
        y = "Snowmelt (mm/h)",
        fill = NULL
      )
  })

  output$return_curve <- renderPlot({
    b <- cell_bits()
    rl <- runoff_level()
    rp <- exp(seq(log(2), log(500), length.out = 200))
    lev <- mixture_tail_return_level(b$mt, rp * b$nep)
    d <- tibble::tibble(rp = rp, level = lev) |> dplyr::filter(is.finite(level))
    validate(need(nrow(d) > 1, "No usable return curve for this cell."))

    ggplot(d, aes(rp, level)) +
      geom_line(colour = ros_palette$accent, linewidth = 0.9) +
      {
        if (is.finite(rl$z)) {
          geom_point(
            data = tibble::tibble(rp = rl$rp, level = rl$z),
            colour = ros_palette$warm,
            size = 3
          )
        }
      } +
      scale_x_log10() +
      labs(
        title = "Runoff marginal for this cell",
        subtitle = "Closed-form mixture tail; the dot is the slider's return period",
        x = "Return period (years)",
        y = "Runoff (mm/h)"
      )
  })

  output$cell_meta <- renderText({
    b <- cell_bits()
    rl <- runoff_level()
    paste(
      c(
        sprintf("cell_id        %s", input$cell_id),
        sprintf("events / year  %.2f", b$nep),
        sprintf("tail components %d", length(b$mt$threshold)),
        sprintf("valid above    %.3f mm/h", mixture_tail_lower(b$mt)),
        sprintf("T = %.1f y", rl$rp),
        sprintf("runoff level   %s", if (is.finite(rl$z)) sprintf("%.4g mm/h", rl$z) else "below tail region"),
        sprintf("peak hours     %d", nrow(b$peaks))
      ),
      collapse = "\n"
    )
  })
}

shinyApp(ui, server)
