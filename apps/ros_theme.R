# Shared look and feel for the rainonsnow Shiny apps.
#
# Source this near the top of an app, after the library() calls:
#
#   source(file.path(repo_root, "apps", "ros_theme.R"))
#   ...
#   ui <- fluidPage(theme = ros_bs_theme(), ...)
#
# Everything degrades gracefully: if bslib is not installed the apps keep their
# default Bootstrap styling rather than failing to start, and the ggplot theme
# is applied through theme_set() so individual plots need no changes.

# Palette. Cool blues for runoff/water, warm tones for the fitted or modelled
# quantity, so the same distinction reads the same way in every app.
ros_palette <- list(
  ink = "#1b2733",
  muted = "#5a6b7b",
  line = "#d7dee5",
  panel = "#ffffff",
  page = "#f6f8fa",
  accent = "#1f6f8b",
  accent_soft = "#8fc7d8",
  warm = "#c1663a",
  warm_soft = "#e8b18f",
  good = "#3f7d5c",
  bad = "#a63d40"
)

# Discrete colours for model comparisons (forest vs GP tail vs naive, etc.).
ros_model_colours <- c(
  "Random Forest" = ros_palette$accent,
  "GP conversion" = ros_palette$warm,
  "Forest mixture" = ros_palette$accent,
  "Forest mixture + GP tail" = ros_palette$warm,
  "Empirical" = ros_palette$muted,
  "Naive POT" = ros_palette$good,
  "Copula transport" = ros_palette$good
)

#' Bootstrap theme for the apps. Returns NULL when bslib is unavailable, which
#' `fluidPage(theme = NULL)` accepts.
ros_bs_theme <- function() {
  if (!requireNamespace("bslib", quietly = TRUE)) {
    return(NULL)
  }
  bslib::bs_theme(
    version = 5,
    bg = ros_palette$page,
    fg = ros_palette$ink,
    primary = ros_palette$accent,
    secondary = ros_palette$muted,
    success = ros_palette$good,
    danger = ros_palette$bad,
    # Stick to stacks that are present on essentially any machine, so the apps
    # look the same offline as they do with web fonts available.
    base_font = bslib::font_collection(
      bslib::font_google("Inter", local = FALSE),
      "Segoe UI",
      "Helvetica Neue",
      "Arial",
      "sans-serif"
    ),
    heading_font = bslib::font_collection(
      bslib::font_google("Inter", local = FALSE),
      "Segoe UI",
      "Helvetica Neue",
      "Arial",
      "sans-serif"
    ),
    code_font = c("SFMono-Regular", "Menlo", "Consolas", "monospace"),
    "body-color" = ros_palette$ink,
    "border-color" = ros_palette$line
  )
}

#' A little CSS the Bootstrap variables do not cover: tighter sidebars, calmer
#' plot borders, and a subtitle style for the one-line explanations that sit
#' under each app title.
ros_css <- function() {
  shiny::tags$style(shiny::HTML(sprintf(
    "
    body { background: %s; }
    .well, .card, .panel {
      background: %s;
      border: 1px solid %s;
      border-radius: 10px;
      box-shadow: 0 1px 2px rgba(16,24,32,.04);
    }
    h1, h2, h3, h4 { letter-spacing: -0.01em; font-weight: 600; }
    h2.ros-title { font-size: 1.5rem; margin: 0 0 .15rem 0; }
    p.ros-subtitle { color: %s; margin: 0 0 1rem 0; font-size: .92rem; }
    .ros-note {
      color: %s; font-size: .85rem; line-height: 1.45;
      border-left: 3px solid %s; padding: .35rem 0 .35rem .6rem; margin: .6rem 0;
    }
    .shiny-input-container { margin-bottom: .7rem; }
    .form-control, .selectize-input { border-color: %s; }
    .shiny-plot-output { background: %s; border-radius: 8px; }
    .nav-tabs .nav-link.active { font-weight: 600; }
    ",
    ros_palette$page,
    ros_palette$panel,
    ros_palette$line,
    ros_palette$muted,
    ros_palette$muted,
    ros_palette$accent_soft,
    ros_palette$line,
    ros_palette$panel
  )))
}

#' Standard header: title plus a one-line description of what the app is for.
ros_header <- function(title, subtitle = NULL) {
  shiny::tagList(
    ros_css(),
    shiny::h2(title, class = "ros-title"),
    if (!is.null(subtitle)) shiny::p(subtitle, class = "ros-subtitle")
  )
}

#' A short explanatory note, for stating what a panel is actually showing.
ros_note <- function(...) {
  shiny::div(class = "ros-note", ...)
}

#' Consistent ggplot theme. Lighter than theme_bw(), with the panel border kept
#' because most of these plots are read quantitatively.
ros_ggtheme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = ros_palette$ink),
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.05),
      plot.subtitle = ggplot2::element_text(
        colour = ros_palette$muted,
        size = base_size * 0.9
      ),
      plot.caption = ggplot2::element_text(
        colour = ros_palette$muted,
        size = base_size * 0.78
      ),
      axis.title = ggplot2::element_text(colour = ros_palette$muted),
      axis.text = ggplot2::element_text(colour = ros_palette$muted),
      panel.grid.major = ggplot2::element_line(
        colour = ros_palette$line,
        linewidth = 0.3
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = ros_palette$line,
        fill = NA,
        linewidth = 0.5
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        colour = ros_palette$ink,
        size = base_size * 0.9
      ),
      strip.background = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.key.height = ggplot2::unit(0.8, "lines"),
      plot.background = ggplot2::element_rect(
        fill = ros_palette$panel,
        colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = ros_palette$panel,
        colour = NA
      )
    )
}

#' Map theme: no axes, no grid -- the coordinates are not the point.
ros_ggtheme_map <- function(base_size = 12) {
  ros_ggtheme(base_size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )
}

#' Colour scale for named models, falling back to the default for unknown names.
ros_scale_model <- function(...) {
  ggplot2::scale_colour_manual(values = ros_model_colours, na.value = ros_palette$muted, ...)
}

# Apply the theme to every plot in the app that sources this file.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::theme_set(ros_ggtheme())
}
