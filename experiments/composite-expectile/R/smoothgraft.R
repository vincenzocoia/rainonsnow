# ---------------------------------------------------------------------------
# Return levels from a smooth graft: empirical body, fitted GEV tail.
#
# Uses the real construction from probaverse/distplyr (branch
# claude/distplyr-smooth-graft-spond9) -- the hazard-mixture smooth graft --
# rather than a reimplementation. The graft's survival is
#
#   S(x) = M(x) C(x),   M = (1-w) Sbody + w Stail,
#   C(x) = exp{ -int_{x0}^{x} w'(s) [Stail(s) - Sbody(s)] / M(s) ds }.
#
# With a weight that is exactly 0 below xlo and exactly 1 above xhi, the
# construction collapses to three regimes, which is what makes it cheap to
# evaluate on a return-period grid:
#
#   x <  xlo : w = 0, w' = 0, so C = 1 and S = Sbody. The graft IS the
#              empirical distribution.
#   xlo..xhi : the transition; needs the integral.
#   x >  xhi : w = 1, w' = 0, so C is frozen at c = C(xhi) and
#              S(x) = c * Stail(x). The graft is the fitted tail, rescaled by a
#              single constant that reconciles it with the body's mass.
#
# Both limits and the constancy of c are verified in scripts/98-validate.R
# against distplyr's own eval_survival and eval_quantile.
# ---------------------------------------------------------------------------

# The composite weight transported to the x-axis: the same smoothstep, running
# over the interval between the sample's p0- and p1-quantiles, so the handover
# occupies the same range of levels that the fitting weight does.
graft_weight <- function(y, p0, p1) {
  xlo <- as.numeric(stats::quantile(y, p0, type = 1))
  xhi <- as.numeric(stats::quantile(y, p1, type = 1))
  if (xhi <= xlo) xhi <- xlo + 1e-8 * max(1, abs(xlo))
  list(
    xlo = xlo, xhi = xhi,
    w = function(x) { s <- pmin(1, pmax(0, (x - xlo) / (xhi - xlo))); s^2 * (3 - 2 * s) },
    wp = function(x) { s <- pmin(1, pmax(0, (x - xlo) / (xhi - xlo)))
                       ifelse(s <= 0 | s >= 1, 0, 6 * s * (1 - s) / (xhi - xlo)) }
  )
}

# Return levels of the smooth graft at exceedance probabilities `ex` (= 1/T).
graft_return_levels <- function(y, theta, p0, p1, ex) {
  if (any(is.na(theta))) return(rep(NA_real_, length(ex)))
  gw <- graft_weight(y, p0, p1)
  body <- distionary::dst_empirical(y)
  tail <- distionary::dst_gev(theta[1], theta[2], theta[3])
  g <- try(distplyr::smooth_graft_right(body, tail, weight = gw$w,
                                        weight_deriv = gw$wp), silent = TRUE)
  if (inherits(g, "try-error")) return(rep(NA_real_, length(ex)))

  # c, read off exactly at the top of the handover, where w = 1 so that
  # M(xhi) = Stail(xhi) and hence S(xhi) = c * Stail(xhi). Evaluating here
  # rather than further out matters: a fit with a negative shape has a finite
  # upper endpoint, and a probe beyond it returns Stail = 0 and no usable c.
  s_tail_hi <- gev_survival(gw$xhi, theta[1], theta[2], theta[3])
  if (!is.finite(s_tail_hi) || s_tail_hi <= 0) {
    # The fitted tail cannot even reach the top of the handover: its upper
    # endpoint lies at or below the sample's p1-quantile. The construction has
    # no tail to hand over to, so there is no graft for this replicate.
    return(structure(rep(NA_real_, length(ex)), reason = "tail endpoint below handover"))
  }
  cst <- as.numeric(distionary::eval_survival(g, at = gw$xhi)) / s_tail_hi
  if (!is.finite(cst) || cst <= 0) return(rep(NA_real_, length(ex)))

  e_body <- as.numeric(distionary::eval_survival(body, at = gw$xlo))  # start of handover
  e_tail <- cst * s_tail_hi                                           # end of handover

  out <- numeric(length(ex))
  for (i in seq_along(ex)) {
    e <- ex[i]
    if (e >= e_body) {
      # Body regime: the graft is the empirical distribution exactly.
      out[i] <- as.numeric(distionary::eval_quantile(body, at = 1 - e))
    } else if (e <= e_tail) {
      # Tail regime: S = c * Stail, so invert the fitted GEV at the rescaled level.
      out[i] <- qgev(1 - e / cst, theta[1], theta[2], theta[3])
    } else {
      # Transition: the only regime that needs the correction integral.
      f <- function(x) as.numeric(distionary::eval_survival(g, at = x)) - e
      out[i] <- tryCatch(stats::uniroot(f, c(gw$xlo, gw$xhi), tol = 1e-9)$root,
                         error = function(...) NA_real_)
    }
  }
  attr(out, "c") <- cst
  out
}
