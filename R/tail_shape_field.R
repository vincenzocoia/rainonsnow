# Sharing tail-shape information between grid cells.
#
# fit_gpd_shared_shape() removes the scatter *within* a cell by giving all of
# that cell's peak hours one shape. The per-cell shape estimates that come out
# are still noisy between cells, and the shape of a rain-on-snow runoff tail is
# a smooth physical property of the terrain rather than something that should
# jump between adjacent 0.1-degree cells. So the natural next step is to borrow
# strength across neighbours.

#' Neighbours of each cell on a regular grid
#'
#' Cells are neighbours when they are within `radius` grid steps in both
#' directions (a Chebyshev ball, so `radius = 1` is the 8 surrounding cells).
#' The step is inferred from the smallest positive spacing on each axis.
#'
#' @param x,y Numeric vectors of cell centre coordinates.
#' @param radius Neighbourhood radius in grid steps.
#' @param include_self Whether a cell is its own neighbour.
#' @returns A list of integer vectors of indices, one per cell.
#' @examples
#' g <- expand.grid(x = 1:3, y = 1:3)
#' grid_neighbours(g$x, g$y)[[5]]
#' @export
grid_neighbours <- function(x, y, radius = 1, include_self = FALSE) {
  step <- function(v) {
    u <- sort(unique(v))
    if (length(u) < 2L) {
      return(1)
    }
    d <- diff(u)
    min(d[d > 0])
  }
  dx <- step(x)
  dy <- step(y)
  ix <- round(x / dx)
  iy <- round(y / dy)

  n <- length(x)
  # Hash grid coordinates so lookup is O(n * (2r+1)^2) rather than O(n^2).
  key <- paste(ix, iy, sep = ",")
  idx <- stats::setNames(seq_len(n), key)
  offsets <- expand.grid(
    dx = seq(-radius, radius),
    dy = seq(-radius, radius)
  )
  if (!include_self) {
    offsets <- offsets[!(offsets$dx == 0 & offsets$dy == 0), , drop = FALSE]
  }

  lapply(seq_len(n), function(i) {
    k <- paste(ix[i] + offsets$dx, iy[i] + offsets$dy, sep = ",")
    out <- idx[k]
    unname(out[!is.na(out)])
  })
}

#' Smooth per-cell tail shapes toward their neighbours
#'
#' Empirical-Bayes shrinkage of each cell's shape estimate toward the mean of
#' its neighbours. A cell whose own estimate is precise moves very little; a
#' cell with a noisy estimate is pulled most of the way to its neighbourhood.
#'
#' Writing `xi_c` for the estimate in cell `c`, `s_c` for its standard error and
#' `m_c` for the neighbourhood mean, the smoothed value is the precision-weighted
#' combination
#' \deqn{(xi_c / s_c^2 + m_c / tau^2) / (1 / s_c^2 + 1 / tau^2).}
#' The between-cell variance `tau^2` is the DerSimonian--Laird estimator, the
#' standard choice for this problem: the spread left over once sampling noise is
#' removed, floored at zero. The naive `var(xi) - mean(s^2)` is not used because
#' `mean(s^2)` is dominated by the least precise cells, which drags `tau^2` down
#' and over-smooths everything else. When `tau^2` is zero the shapes are
#' statistically indistinguishable and every cell collapses to the common value,
#' which is the right answer in that case.
#'
#' @param shape Numeric vector of per-cell shape estimates.
#' @param se Numeric vector of standard errors, e.g. `shape_se` from
#'   [fit_gpd_shared_shape()].
#' @param x,y Cell coordinates, passed to [grid_neighbours()].
#' @param radius Neighbourhood radius in grid steps.
#' @param min_neighbours Cells with fewer usable neighbours than this fall back
#'   to the global mean.
#' @returns A list with `shape` (the smoothed estimates), `tau` (the estimated
#'   between-cell standard deviation), and `weight` (how much of each cell's own
#'   estimate was retained, between 0 and 1).
#' @examples
#' set.seed(1)
#' g <- expand.grid(x = 1:5, y = 1:5)
#' truth <- 0.1 + 0.02 * g$x
#' est <- truth + stats::rnorm(nrow(g), sd = 0.15)
#' out <- smooth_tail_shape(est, se = rep(0.15, nrow(g)), x = g$x, y = g$y)
#' # Smoothing moves the estimates back toward the underlying field.
#' c(raw = mean((est - truth)^2), smoothed = mean((out$shape - truth)^2))
#' @export
smooth_tail_shape <- function(
  shape,
  se,
  x,
  y,
  radius = 1,
  min_neighbours = 2L
) {
  n <- length(shape)
  se <- rep_len(se, n)
  usable <- is.finite(shape) & is.finite(se) & se > 0

  if (sum(usable) < 2L) {
    return(list(shape = shape, tau = NA_real_, weight = rep(1, n)))
  }

  nb <- grid_neighbours(x, y, radius = radius, include_self = FALSE)
  # Neighbours are combined by precision, not equally: a neighbour whose own
  # shape is barely determined should not pull as hard as a well-determined one.
  prec <- ifelse(usable, 1 / se^2, 0)
  global <- sum(shape[usable] * prec[usable]) / sum(prec[usable])

  neighbourhood_mean <- function(i) {
    j <- nb[[i]]
    j <- j[usable[j]]
    if (length(j) >= min_neighbours) {
      sum(shape[j] * prec[j]) / sum(prec[j])
    } else {
      global
    }
  }
  prior <- vapply(seq_len(n), function(i) {
    if (usable[i]) neighbourhood_mean(i) else NA_real_
  }, numeric(1L))

  # The tau^2 = 0 limit means the cell and its neighbours share one shape, so
  # the best estimate pools them ALL -- the cell's own data included. Using the
  # neighbours-only mean there would throw that cell's observations away.
  nb_self <- grid_neighbours(x, y, radius = radius, include_self = TRUE)
  prior_pooled <- vapply(seq_len(n), function(i) {
    if (!usable[i]) {
      return(NA_real_)
    }
    j <- nb_self[[i]]
    j <- j[usable[j]]
    if (length(j) > min_neighbours) {
      sum(shape[j] * prec[j]) / sum(prec[j])
    } else {
      global
    }
  }, numeric(1L))

  # Between-cell variance, DerSimonian--Laird, measured around the LOCAL prior
  # rather than the global mean. The shrinkage target is the neighbourhood, so
  # the relevant dispersion is how far a cell sits from its neighbours; on a
  # field with any spatial trend the global spread is much larger than that and
  # would leave the estimates barely smoothed at all.
  yi <- shape[usable]
  wi <- 1 / se[usable]^2
  resid <- yi - prior[usable]
  s1 <- sum(wi)
  s2 <- sum(wi^2)
  q <- sum(wi * resid^2)
  denom <- s1 - s2 / s1
  tau2 <- if (denom > 0) max(0, (q - (length(yi) - 1)) / denom) else 0
  tau <- sqrt(tau2)

  out <- shape
  wt <- rep(1, n)
  for (i in seq_len(n)) {
    if (!usable[i]) {
      next
    }
    if (tau2 <= 0) {
      # No detectable variation between cells beyond sampling noise, so the
      # cells are treated as sharing one shape.
      out[i] <- prior_pooled[i]
      wt[i] <- 0
      next
    }
    w <- (1 / se[i]^2) / (1 / se[i]^2 + 1 / tau2)
    out[i] <- w * shape[i] + (1 - w) * prior[i]
    wt[i] <- w
  }
  list(shape = out, tau = tau, weight = wt)
}
