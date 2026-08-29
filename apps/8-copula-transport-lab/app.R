# Copula transport lab: recovering a marginal tail from conditional tails.
#
# Run from the package root:
#   shiny::runApp("apps/8-copula-transport-lab")
#
# Needs no derived data. Everything on screen is simulated from a
# data-generating process you choose in the sidebar, so the TRUTH is known
# exactly and every estimate can be scored against it.
#
# WHAT IT IS FOR
#
# The marginal of Y cannot be recovered by mixing the learned conditionals over
# the covariate values in the sample: they all share the conditional tail index,
# so the mixture is regularly varying with THAT index, lighter than the
# marginal's by the factor CEVI. The copula supplies the missing step. With
# h(v|u) = dC(u,v)/du,
#
#     S_Y(y) = h^{-1}( S_{Y|X=x}(y) | u ),      u = F_X(x),
#
# so every single covariate value gives an estimate of the whole marginal. This
# app builds those estimates, combines them, and shows the answer against the
# truth and against the obvious alternative -- a generalized Pareto fitted to Y
# on its own, ignoring X.
#
# WHY THE DGP IS BUILT THE WAY IT IS
#
# F_X and Y | X are chosen; the marginal of Y is DERIVED from them. Nothing
# about the marginal is picked for convenience, and it lands in no canned
# family: the "DGP" tab shows its local tail index still moving at exceedance
# probabilities far below anything a sample reaches. So the canned fit is not
# imprecise, it is fitting the wrong shape of curve, and more data does not
# rescue it.
#
# Suggests: shiny, ggplot2

library(shiny)
library(ggplot2)

repo_root <- here::here()
devtools::load_all(repo_root, quiet = TRUE)
source(file.path(repo_root, "apps", "ros_theme.R"))
theme_set(ros_ggtheme())

PROBS <- c(1e-2, 1e-3, 1e-4, 1e-5)

# ---------------------------------------------------------------------------
ui <- fluidPage(
  theme = ros_bs_theme(),
  ros_header(
    "Copula transport lab",
    paste(
      "Every covariate value transports its own conditional tail back into an",
      "estimate of the whole marginal. Combine them, and score against a truth",
      "that is derived rather than assumed."
    )
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h5("The process"),
      selectInput(
        "x_dist",
        "Covariate distribution",
        c(
          "Pareto  (heavier than Y|X, so CEVI < 1)" = "pareto",
          "LogNormal  (all moments, so CEVI = 1)" = "lognormal"
        )
      ),
      conditionalPanel(
        "input.x_dist == 'pareto'",
        sliderInput("alpha", "Pareto index of X", 1.5, 12, 4, step = 0.25)
      ),
      conditionalPanel(
        "input.x_dist == 'lognormal'",
        sliderInput("sdlog", "Log-scale spread of X", 0.25, 3, 1, step = 0.05)
      ),
      sliderInput("xi_cond", "Shape of Y | X", 0.01, 0.6, 0.1, step = 0.01),
      sliderInput("n", "Sample size", 250, 20000, 3000, step = 250),
      actionButton("reseed", "New sample", class = "btn-sm btn-primary"),
      hr(),
      h5("Learning the conditionals"),
      sliderInput("bandwidth", "Kernel bandwidth (rank scale)", 0.01, 0.2,
                  0.03, step = 0.005),
      sliderInput("n_anchor", "Covariate values", 5, 80, 30, step = 5),
      sliderInput("tail_frac", "Upper fraction used for the tail", 0.05, 0.7,
                  0.35, step = 0.05),
      checkboxInput("shared", "One tail shape across covariate values", TRUE),
      hr(),
      h5("Transport and combine"),
      selectInput(
        "copula",
        "Copula used for the transport",
        c(
          "the true one, from the DGP" = "true",
          "Gaussian, fitted" = "gaussian",
          "survival Clayton, fitted" = "survival_clayton",
          "Gaussian, set by hand" = "gaussian_manual"
        )
      ),
      conditionalPanel(
        "input.copula == 'gaussian_manual'",
        sliderInput("rho", "Correlation", -0.95, 0.95, 0.7, step = 0.05)
      ),
      radioButtons(
        "combine",
        "Combination rule",
        c("pointwise median" = "median",
          "mean" = "mean",
          "weighted mean" = "weighted_mean")
      )
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Marginal tail",
          br(),
          plotOutput("survival_plot", height = "460px"),
          ros_note(
            "One thin line per covariate value: each is that covariate value's",
            "own estimate of the whole marginal, obtained by inverting the",
            "copula's h-function on its conditional tail. They all estimate the",
            "SAME curve, so the spread between them is a free diagnostic of how",
            "far the conditional model and the copula are from agreeing."
          ),
          h5("Return levels"),
          tableOutput("levels")
        ),
        tabPanel(
          "Across the covariate",
          br(),
          plotOutput("by_x_plot", height = "520px"),
          ros_note(
            "Top: each covariate value's estimate of one return level, against",
            "the truth. Bottom: how deep an extrapolation that estimate needs --",
            "the number of log units between the conditional's own threshold and",
            "the conditional exceedance probability the target event corresponds",
            "to. The estimate is y proportional to s^-xi across that span, so a",
            "small error in the shape is a large error in the answer wherever",
            "the span is wide. That is what CEVI < 1 buys: it narrows the span."
          )
        ),
        tabPanel(
          "The process",
          br(),
          plotOutput("dgp_plot", height = "480px"),
          ros_note(
            "The marginal's local tail index, -dlog S / dlog y, against the",
            "asymptote Breiman's lemma gives it. Where the two differ, a",
            "generalized Pareto fitted to Y alone reads the local slope and then",
            "extrapolates it as if it had settled. The conditional has no such",
            "problem: it is exactly generalized Pareto everywhere, by",
            "construction."
          ),
          verbatimTextOutput("dgp_text")
        ),
        tabPanel(
          "What this is",
          br(),
          uiOutput("about")
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  dgp <- reactive({
    transport_dgp(
      x_dist = input$x_dist,
      sdlog = input$sdlog %||% 1,
      alpha = input$alpha %||% 4,
      xi_cond = input$xi_cond
    )
  })

  sample_data <- reactive({
    input$reseed
    d <- dgp()
    isolate(td_simulate(d, input$n))
  })

  conditionals <- reactive({
    s <- sample_data()
    td_conditionals(
      s$u,
      s$y,
      n_anchor = input$n_anchor,
      bandwidth = input$bandwidth,
      tail_frac = input$tail_frac
    )
  })

  # The copula the transport is run through. "true" is the exact h-inverse the
  # DGP implies, supplied as a function; the others are fitted from the sample's
  # Kendall tau, so the difference between them is model risk in isolation.
  transport_family <- reactive({
    s <- sample_data()
    d <- dgp()
    tau <- stats::cor(s$u, s$y, method = "kendall")
    switch(
      input$copula,
      true = list(family = function(sc, u) td_transport(d, sc, u), par = NA,
                  label = "true copula"),
      gaussian = list(family = "gaussian", par = sin(pi / 2 * tau),
                      label = sprintf("Gaussian, rho = %.2f",
                                      sin(pi / 2 * tau))),
      survival_clayton = list(
        family = "survival_clayton",
        par = max(2 * tau / (1 - tau), 0.05),
        label = sprintf("survival Clayton, theta = %.2f",
                        max(2 * tau / (1 - tau), 0.05))
      ),
      gaussian_manual = list(family = "gaussian", par = input$rho,
                             label = sprintf("Gaussian, rho = %.2f", input$rho))
    )
  })

  ensemble <- reactive({
    cc <- conditionals()
    tf <- transport_family()
    validate(need(
      any(!vapply(cc$pieces, is.null, logical(1L))),
      "No conditional has enough exceedances. Widen the bandwidth."
    ))
    transport_ensemble(
      cc$pieces,
      cc$u,
      family = tf$family,
      par = if (is.na(tf$par)) 0.5 else tf$par,
      shared_shape = input$shared
    )
  })

  canned <- reactive(td_canned_gp(sample_data()$y))

  # --- the marginal tail ---------------------------------------------------
  output$survival_plot <- renderPlot({
    d <- dgp()
    me <- ensemble()
    g <- canned()
    lo <- marginal_ensemble_lower(me)
    hi <- td_marginal_quantile(d, 1e-6)
    yy <- exp(seq(log(lo), log(max(hi, lo * 10)), length.out = 250L))

    paths <- marginal_ensemble_paths(me, yy)
    path_df <- data.frame(
      y = rep(yy, ncol(paths)),
      s = as.vector(paths),
      u = rep(round(me$u, 3), each = length(yy))
    )
    combined <- data.frame(
      y = yy,
      s = marginal_ensemble_survival(me, yy, input$combine),
      what = "transport, combined"
    )
    truth <- data.frame(y = yy, s = td_marginal_survival(d, yy), what = "truth")
    cannedd <- data.frame(y = yy, s = g$survival(yy),
                          what = "generalized Pareto on Y alone")

    s <- sample_data()
    ys <- sort(s$y[s$y > lo], decreasing = TRUE)
    emp <- data.frame(
      y = ys,
      s = seq_along(ys) / length(s$y)
    )

    ggplot(mapping = aes(y, s)) +
      geom_line(
        data = path_df,
        aes(group = u),
        colour = ros_palette$accent_soft,
        alpha = 0.55,
        linewidth = 0.35
      ) +
      geom_point(data = emp, colour = ros_palette$muted, size = 0.7,
                 alpha = 0.5) +
      geom_line(
        data = rbind(combined, truth, cannedd),
        aes(colour = what, linetype = what),
        linewidth = 0.9
      ) +
      scale_x_log10() +
      scale_y_log10(labels = function(z) format(z, scientific = TRUE)) +
      scale_colour_manual(values = c(
        "truth" = ros_palette$ink,
        "transport, combined" = ros_palette$accent,
        "generalized Pareto on Y alone" = ros_palette$warm
      )) +
      scale_linetype_manual(values = c(
        "truth" = "dashed",
        "transport, combined" = "solid",
        "generalized Pareto on Y alone" = "dotted"
      )) +
      labs(
        x = "y", y = "P(Y > y)", colour = NULL, linetype = NULL,
        title = "Marginal survival, one estimate per covariate value",
        subtitle = sprintf(
          "%s | %s | conditional shape %s",
          transport_family()$label,
          paste0(input$combine, " across ", length(me$u), " covariate values"),
          if (input$shared) {
            sprintf("%.3f shared (true %.3f)", me$shape[1], d$xi_cond)
          } else {
            sprintf("%.3f to %.3f (true %.3f)", min(me$shape), max(me$shape),
                    d$xi_cond)
          }
        )
      ) +
      theme(legend.position = "top")
  })

  output$levels <- renderTable(
    {
      d <- dgp()
      me <- ensemble()
      g <- canned()
      truth <- td_marginal_quantile(d, PROBS)
      tr <- marginal_ensemble_quantile(me, PROBS, input$combine)
      cn <- g$threshold +
        gpd_quantile_upper(PROBS / g$tail_prob, g$scale, g$shape)
      data.frame(
        `exceedance probability` = format(PROBS, scientific = TRUE),
        truth = truth,
        transport = tr,
        `transport / truth` = tr / truth,
        `GP on Y alone` = cn,
        `GP / truth` = cn / truth,
        check.names = FALSE
      )
    },
    digits = 3,
    width = "100%"
  )

  # --- across the covariate ------------------------------------------------
  output$by_x_plot <- renderPlot({
    d <- dgp()
    me <- ensemble()
    target <- 1e-4
    truth <- td_marginal_quantile(d, target)

    per_x <- vapply(
      seq_along(me$u),
      function(j) {
        one <- me
        for (f in c("u", "threshold", "tail_prob", "scale", "shape",
                    "weights")) {
          one[[f]] <- me[[f]][j]
        }
        marginal_ensemble_quantile(one, target)
      },
      numeric(1L)
    )
    # How far past its own threshold each conditional has to be pushed for the
    # target event, in log units.
    s_cond <- me$tail_prob *
      gpd_survival(truth - me$threshold, me$scale, me$shape)
    span <- log10(me$tail_prob / pmax(s_cond, 1e-300))

    top <- ggplot(data.frame(u = me$u, q = per_x), aes(u, q)) +
      geom_hline(yintercept = truth, linetype = "dashed",
                 colour = ros_palette$ink) +
      geom_line(colour = ros_palette$accent, linewidth = 0.7) +
      geom_point(colour = ros_palette$accent, size = 1.6) +
      geom_hline(
        yintercept = marginal_ensemble_quantile(me, target, input$combine),
        colour = ros_palette$warm
      ) +
      labs(
        x = NULL, y = "estimated 1e-4 level",
        title = "Each covariate value's estimate of the same return level",
        subtitle = paste(
          "dashed: the truth.  orange: the combined estimate.",
          "Every point is an estimate of the same number."
        )
      )

    bottom <- ggplot(data.frame(u = me$u, span = span), aes(u, span)) +
      geom_col(fill = ros_palette$accent_soft) +
      labs(
        x = "u = F_X(x)",
        y = "log10 units extrapolated",
        subtitle = paste(
          "How far past its own threshold each conditional is pushed to reach",
          "the target. Narrower is better; this is what CEVI < 1 buys."
        )
      )

    if (requireNamespace("patchwork", quietly = TRUE)) {
      patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(2, 1))
    } else {
      gridExtra_stack(top, bottom)
    }
  })

  # --- the process ---------------------------------------------------------
  output$dgp_plot <- renderPlot({
    d <- dgp()
    ps <- 10^seq(-0.3, -12, length.out = 120L)
    yy <- td_marginal_quantile(d, ps)
    df <- data.frame(
      s = ps,
      evi = td_local_evi(d, yy),
      what = "marginal of Y (derived)"
    )
    ggplot(df, aes(s, evi)) +
      geom_line(colour = ros_palette$accent, linewidth = 0.9) +
      geom_hline(
        yintercept = td_asymptotic_evi(d),
        linetype = "dashed",
        colour = ros_palette$ink
      ) +
      geom_hline(
        yintercept = d$xi_cond,
        linetype = "dotted",
        colour = ros_palette$warm,
        linewidth = 0.8
      ) +
      annotate(
        "rect",
        xmin = 1 / input$n, xmax = 0.5, ymin = -Inf, ymax = Inf,
        fill = ros_palette$muted, alpha = 0.08
      ) +
      scale_x_log10(labels = function(z) format(z, scientific = TRUE)) +
      labs(
        x = "P(Y > y)",
        y = "local tail index",
        title = "The marginal's local tail index against its asymptote",
        subtitle = paste0(
          "dashed: the asymptotic marginal index.  dotted: the conditional's, ",
          "which is exact at every y.  Shaded: what a sample of ",
          input$n, " can see."
        )
      )
  })

  output$dgp_text <- renderText({
    d <- dgp()
    out <- utils::capture.output(print(d))
    ps <- 10^-(1:6)
    yy <- td_marginal_quantile(d, ps)
    c(
      out,
      "",
      sprintf("%12s %12s %12s", "S_Y(y)", "y", "local EVI"),
      sprintf("%12.0e %12.3f %12.3f", ps, yy, td_local_evi(d, yy))
    ) |>
      paste(collapse = "\n")
  })

  output$about <- renderUI({
    HTML(
      "<p>The marginal of <code>Y</code> cannot be recovered by mixing the
      learned conditionals over the covariate values in the sample. They all
      share the conditional tail index, so any finite mixture of them is
      regularly varying with <em>that</em> index, lighter than the marginal's by
      the factor CEVI. The marginal only gets its full weight from covariate
      values beyond any sample, and no amount of learning reaches them.</p>

      <p>The copula does. Writing <code>h(v|u) = dC(u,v)/du</code>, the
      conditional and the marginal are related by
      <code>S_{Y|X=x}(y) = h(S_Y(y) | u)</code>, so inverting <code>h</code> in
      its first argument turns any <em>one</em> conditional into an estimate of
      the whole marginal. There is one such estimate per covariate value.</p>

      <p><b>Combining them is an averaging problem, not a mixing one.</b> That
      distinction decides which operations are legal. Mixing conditionals is the
      law of total probability and admits only the weighted mean. These are
      repeated estimates of one quantity, so the pointwise <em>median</em> is
      available -- and it is well defined as a distribution, since a pointwise
      median of monotone survival curves is itself a monotone survival curve.
      Turn off the shared shape in the sidebar to see why it matters: the mean
      is dragged upward by whichever covariate value drew the heaviest shape,
      which is the same domination that makes an unpooled mixture inherit
      <code>max(xi_i)</code>, reappearing one level up.</p>

      <p><b>The precondition is CEVI &lt; 1.</b> Switch the covariate to
      LogNormal and the marginal inherits the conditional's own tail index, so
      conditioning buys no reduction in extremity: reaching a 1e-4 marginal
      event from the middle of the covariate range still needs a conditional
      event around 1e-8, and the &quot;log units extrapolated&quot; panel shows
      the span the estimate has to cross. Since the estimate is
      <code>y</code> proportional to <code>s^-xi</code> over that span, a shape
      error of 0.01 across eighteen log units is a factor of 1.2 in the answer.
      Where conditioning genuinely makes an extreme less extreme the span is
      short and the transport is worth doing; where it does not, there is
      nothing to decompose, and the ensemble spread says so before you are
      misled.</p>

      <p>The truth on every plot is derived from the chosen <code>F_X</code> and
      <code>Y|X</code> by quadrature, never fitted, and the marginal it produces
      belongs to no canned family &mdash; which is what makes the
      &quot;generalized Pareto on Y alone&quot; comparison a fair fight rather
      than a rigged one.</p>"
    )
  })
}

# Fallback stacker when patchwork is not installed.
gridExtra_stack <- function(a, b) {
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(a, b, ncol = 1, heights = c(2, 1))
  } else {
    a
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x

shinyApp(ui, server)
