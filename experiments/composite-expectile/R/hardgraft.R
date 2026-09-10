# ---------------------------------------------------------------------------
# Hard graft: empirical body below a threshold, parametric tail above.
#
# The tail is rescaled so the survival function is continuous at the join,
#     S(x) = (1 - v) S_theta(x) / S_theta(u)   for x > u,   u = Fhat^-1(v),
# which is the ordinary peaks-over-threshold construction. The cdf is
# continuous but the density is not: there is a kink at u, and the fitted tail
# shape is imposed abruptly rather than blended in.
#
# This is the comparison the smooth graft is meant to beat. If a hard graft
# already repairs the body, the smooth handover is machinery without a purpose.
# ---------------------------------------------------------------------------

hard_graft_return_levels <- function(y, theta, ex, family = c("gev", "gpd"), v = 0.90) {
  family <- match.arg(family)
  if (anyNA(theta)) return(rep(NA_real_, length(ex)))
  u  <- as.numeric(stats::quantile(y, v, type = 7))
  Su <- if (family == "gpd") sgpd(u, theta[1], theta[2], theta[3])
        else                 gev_survival(u, theta[1], theta[2], theta[3])
  if (!is.finite(Su) || Su <= 0) return(rep(NA_real_, length(ex)))
  out <- numeric(length(ex))
  body <- ex >= (1 - v)
  if (any(body))
    out[body] <- as.numeric(stats::quantile(y, 1 - ex[body], type = 7))
  if (any(!body)) {
    # solve (1-v) S_theta(x)/S_theta(u) = e
    s_target <- ex[!body] * Su / (1 - v)
    out[!body] <- if (family == "gpd")
      qgpd(1 - s_target, theta[1], theta[2], theta[3])
    else
      qgev(1 - s_target, theta[1], theta[2], theta[3])
  }
  out
}
