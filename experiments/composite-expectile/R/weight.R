# ---------------------------------------------------------------------------
# The weight function w(p) of the composite estimators.
#
# A trapezoid: zero below p0, rising over [p0, p1], one over [p1, p2], falling
# over [p2, p3], zero above p3. Setting p2 = p3 = 1 gives the shape described in
# the proposal -- zero on the body, rising to one as p -> 1.
#
# The falling edge matters more than it looks. Above the level of the largest
# observation the empirical quantile and expectile functions have both saturated
# at the sample maximum, so any weight there can only pull the fitted tail down;
# it is a region the data cannot speak to. scripts/04-shape-diagnostic.R shows
# the fitted shape parameter tracking exactly the share of weight mass that sits
# out there.
# ---------------------------------------------------------------------------

smoothstep <- function(s) { s <- pmin(1, pmax(0, s)); s^2 * (3 - 2 * s) }

make_weight <- function(p0, p1, p2 = 1, p3 = 1) {
  force(p0); force(p1); force(p2); force(p3)
  function(p) {
    up <- smoothstep((p - p0) / (p1 - p0))
    dn <- if (p3 > p2) 1 - smoothstep((p - p2) / (p3 - p2)) else as.numeric(p <= p3)
    up * dn
  }
}
