# Gauss-Legendre nodes and weights on [0, 1], via Newton on the Legendre
# polynomial. Used to integrate the composite loss over levels.
gauss_legendre_01 <- function(n) {
  i <- seq_len(n)
  x <- cos(pi * (i - 0.25) / (n + 0.5))       # Chebyshev starting guess
  for (it in 1:100) {
    p0 <- rep(1, n); p1 <- x
    for (k in 2:n) { p2 <- ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                     p0 <- p1; p1 <- p2 }
    dp <- n * (x * p1 - p0) / (x^2 - 1)
    dx <- p1 / dp
    x <- x - dx
    if (max(abs(dx)) < 1e-15) break
  }
  w <- 2 / ((1 - x^2) * dp^2)
  list(x = (x + 1) / 2, w = w / 2)            # map [-1,1] -> [0,1]
}

# Quadrature rule for  \int_{p0}^{1} w(p) g(p) dp.
#
# Two things make the raw integrand awkward. Near p = 1 it behaves like
# (1-p)^{1-2 xi}, whose derivative is unbounded at p = 1; the substitution
# 1 - p = (1 - p0) s^2 removes that. And the integrand has kinks wherever the
# fitted quantile/expectile crosses an observation, so a single high-order rule
# converges slowly (badly so for the pinball loss, whose kinks are in the first
# derivative). A composite rule -- a modest Gauss-Legendre panel repeated over
# subintervals of s -- localises each kink to one panel and converges quickly.
make_level_grid <- function(p0, n_panel = 24, n_gl = 8) {
  gl <- gauss_legendre_01(n_gl)
  edges <- seq(0, 1, length.out = n_panel + 1)
  h <- 1 / n_panel
  s <- as.vector(outer(gl$x * h, edges[-length(edges)], "+"))
  wq <- rep(gl$w * h, n_panel)
  p <- 1 - (1 - p0) * s^2
  jac <- 2 * (1 - p0) * s                     # |dp/ds|
  list(p = p, w_quad = wq * jac, n_nodes = length(p))
}
