# ---------------------------------------------------------------------------
# The weight function w(p) of the composite estimators.
#
# w is zero on the body, rises smoothly over [p0, p1], and is 1 on [p1, 1).
# p0 is placed where the Gaussian contamination starts to die out and p1 where
# it is effectively gone, so the estimators only see levels at which the GEV is
# (nearly) the correct model.
# ---------------------------------------------------------------------------

make_weight <- function(p0, p1) {
  force(p0); force(p1)
  function(p) {
    s <- (p - p0) / (p1 - p0)
    s <- pmin(1, pmax(0, s))
    s^2 * (3 - 2 * s)          # smoothstep: C^1, zero slope at both ends
  }
}
