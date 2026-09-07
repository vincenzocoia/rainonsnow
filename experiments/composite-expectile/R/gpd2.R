# ---------------------------------------------------------------------------
# Two-parameter GPD fitting: the threshold is KNOWN to be zero, so only the
# scale and shape are estimated.
#
# This is the right shape for the correctly-specified test. A three-parameter
# GPD with an unknown threshold is non-regular -- the likelihood is monotone in
# mu up to min(y), so mu_hat is a boundary estimate converging at rate n -- and
# fitting one makes the composite estimators pay for a parameter the natural
# reference does not estimate. With mu fixed at zero the problem is regular, the
# MLE is efficient, and every estimator in the comparison has the same two
# unknowns.
# ---------------------------------------------------------------------------

# --- references -------------------------------------------------------------
gpd2_mle <- function(y) {
  z <- y[y > 0]
  if (length(z) < 5) return(rep(NA_real_, 2))
  o <- try(optim(c(log(mean(z)), 0.1), gpd_exceed_nll, z = z, method = "Nelder-Mead",
                 control = list(maxit = 2000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 2))
  c(exp(o$par[1]), o$par[2])
}

# L-moments with the threshold known: lambda1 = s/(1-xi) and
# lambda2 = s/((1-xi)(2-xi)), so xi = 2 - lambda1/lambda2 and s = lambda1(1-xi).
gpd2_lmom <- function(y) {
  z <- sort(y[y > 0]); m <- length(z)
  if (m < 5) return(rep(NA_real_, 2))
  b0 <- mean(z); b1 <- sum((seq_len(m) - 1) / (m - 1) * z) / m
  l1 <- b0; l2 <- 2 * b1 - b0
  if (!is.finite(l2) || l2 <= 0) return(rep(NA_real_, 2))
  xi <- 2 - l1 / l2
  s <- l1 * (1 - xi)
  if (!is.finite(s) || s <= 0) return(rep(NA_real_, 2))
  c(s, xi)
}

# --- composite fitters, mu fixed at zero ------------------------------------
# Each is the three-parameter fitter with mu dropped from the parameter vector;
# the loss machinery is untouched.
.gpd2_optim <- function(obj, start_s) {
  par0 <- c(log(max(1e-6, start_s)), 0.15)
  o <- try(optim(par0, obj, method = "Nelder-Mead",
                 control = list(maxit = 3000, reltol = 1e-12)), silent = TRUE)
  if (inherits(o, "try-error") || o$value >= BIG) return(rep(NA_real_, 2))
  o <- optim(o$par, obj, method = "Nelder-Mead",
             control = list(maxit = 3000, reltol = 1e-12))
  if (o$value >= BIG) return(rep(NA_real_, 2))
  c(exp(o$par[1]), o$par[2])
}

fit_gpd2_composite <- function(y, grid, w_fun, type) {
  p <- grid$p; wts <- grid$w_quad * w_fun(p)
  keep <- wts > 0; p <- p[keep]; wts <- wts[keep]
  ys <- sort(y); n <- length(ys)
  C1 <- c(0, cumsum(ys)); C2 <- c(0, cumsum(ys^2)); S1 <- C1[n+1]; S2 <- C2[n+1]
  wp <- wts * p; wq <- wts * (1 - p)
  obj <- function(par) {
    sigma <- exp(par[1]); xi <- par[2]
    if (!is.finite(sigma) || xi < XI_LO || xi > XI_HI) return(BIG)
    tv <- if (type == "quantile") qgpd(p, 0, sigma, xi) else gpd_expectile(p, 0, sigma, xi)
    if (any(!is.finite(tv))) return(BIG)
    j <- findInterval(tv, ys); c1 <- C1[j+1]; c2 <- C2[j+1]
    if (type == "quantile") sum(wts * (p * (S1 - n * tv) - (c1 - j * tv)))
    else {
      lo <- c2 - 2*tv*c1 + j*tv^2; hi <- (S2-c2) - 2*tv*(S1-c1) + (n-j)*tv^2
      sum(wq * lo + wp * hi)
    }
  }
  .gpd2_optim(obj, sd(y))
}

fit_gpd2_elastile <- function(y, grid, w_fun, alpha) {
  p <- grid$p; wts <- grid$w_quad * w_fun(p)
  keep <- wts > 0; p <- p[keep]; wts <- wts[keep]
  ys <- sort(y); n <- length(ys)
  C1 <- c(0, cumsum(ys)); C2 <- c(0, cumsum(ys^2)); S1 <- C1[n+1]; S2 <- C2[n+1]
  st <- gpd2_lmom(y); if (anyNA(st)) st <- c(sd(y), 0.15)
  piece <- function(tv) {
    j <- findInterval(tv, ys); c1 <- C1[j+1]; c2 <- C2[j+1]
    lo <- c2 - 2*tv*c1 + j*tv^2; hi <- (S2-c2) - 2*tv*(S1-c1) + (n-j)*tv^2
    c(L2 = sum(wts*((1-p)*lo + p*hi)), L1 = sum(wts*(p*(S1-n*tv) - (c1-j*tv))))
  }
  s <- piece(gpd_expectile(p, 0, st[1], st[2]))[["L2"]] /
       piece(qgpd(p, 0, st[1], st[2]))[["L1"]]
  if (!is.finite(s) || s <= 0) s <- max(1e-8, sd(y))
  obj <- function(par) {
    sigma <- exp(par[1]); xi <- par[2]
    if (!is.finite(sigma) || xi < XI_LO || xi > XI_HI) return(BIG)
    tv <- gpd_elastile(p, 0, sigma, xi, alpha, s)
    if (any(!is.finite(tv))) return(BIG)
    q <- piece(tv)
    (alpha/s) * q[["L2"]] + (1 - alpha) * q[["L1"]]
  }
  .gpd2_optim(obj, st[1])
}

fit_gpd2_onesided <- function(y, grid, w_fun, k) {
  cc <- k * max(1e-8, IQR(y))
  p <- grid$p; wts <- grid$w_quad * w_fun(p)
  keep <- wts > 0; p <- p[keep]; wts <- wts[keep]
  ys <- sort(y); n <- length(ys)
  C1 <- c(0, cumsum(ys)); C2 <- c(0, cumsum(ys^2)); S1 <- C1[n+1]; S2 <- C2[n+1]
  obj <- function(par) {
    sigma <- exp(par[1]); xi <- par[2]
    if (!is.finite(sigma) || xi <= XI_LO || xi > XI_HI) return(BIG)
    tv <- gpd_onesided(p, 0, sigma, xi, cc)
    if (any(!is.finite(tv))) return(BIG)
    gpd_onesided_loss(tv, p, wts, ys, C1, C2, S1, S2, n, cc)
  }
  .gpd2_optim(obj, sd(y))
}
