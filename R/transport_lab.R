# A data-generating process for experimenting with copula transport.
#
# THE POINT OF THE CONSTRUCTION
#
# The conditional distribution and the covariate distribution are CHOSEN; the
# marginal of Y is whatever they imply. Nothing about the marginal is picked for
# convenience, and it lands outside every canned family:
#
#     X ~ F_X                       (lognormal or Pareto)
#     Y | X = x  ~  GP(scale = x, shape = xi)
#     S_Y(y) = E[ (1 + xi y / X)^(-1/xi) ]        -- derived, by quadrature
#
# So Y is the product X * W with W ~ GP(1, xi) independent of X: a scale
# mixture. Every conditional is EXACTLY generalized Pareto. The marginal is not
# generalized Pareto at any level a data set can reach, and it is not
# generalized extreme value either, so a canned fit to Y on its own is biased
# while the same family fitted to Y | X is exactly right. That is the whole
# argument for decomposing a marginal extreme into its parts, made into
# arithmetic.
#
# HOW BIASED, CONCRETELY
#
# With X ~ LogNormal(0, 1) and xi = 0.2 the marginal's asymptotic tail index is
# also 0.2 -- a lognormal has every moment, so Breiman's lemma gives
# S_Y(y) ~ E[X^{1/xi}] (xi y)^{-1/xi}. But the local slope -dlog S / dlog y gets
# there only around S = 1e-17. Over the range a sample sees it reads:
#
#     S_Y(y)    1e-1   1e-2   1e-3   1e-4   1e-6   1e-10   1e-17
#     local EVI 0.82   0.52   0.40   0.34   0.28   0.22    0.20
#
# A generalized Pareto fitted to the top 10% of n = 5000 returns shape 0.42
# against a truth of 0.20, and its 1-in-10000 level comes out 27% high. The
# multiplier is what makes this happen: the far tail of Y is reached by an
# ordinary W riding an extreme X, and no amount of data pins that down from the
# marginal alone.
#
# THE COPULA IS DERIVED TOO
#
# [td_transport()] is the exact h-inverse of the copula this DGP implies, so the
# transport step can be run with the true dependence and with a fitted
# parametric one, and the difference between them is model risk in isolation.

#' A scale-mixture data-generating process
#'
#' Chooses `F_X` and `Y | X`, and derives everything else. See the file header
#' of `R/transport_lab.R` for why it is built this way.
#'
#' @param x_dist Distribution of the covariate: `"lognormal"` (all moments
#'   finite, so the marginal tail index equals the conditional one and only the
#'   approach to it is slow) or `"pareto"` (heavy enough to raise the marginal
#'   tail index above the conditional one when `alpha < 1 / xi_cond`).
#' @param sdlog Log-scale spread, for `x_dist = "lognormal"`. Larger values push
#'   the marginal further from generalized Pareto over the observable range.
#' @param alpha Pareto tail index of `X`, for `x_dist = "pareto"`.
#' @param xi_cond Shape of the conditional generalized Pareto. Every `Y | X = x`
#'   has this shape exactly.
#' @returns An object of class `transport_dgp`.
#' @examples
#' d <- transport_dgp()
#' td_marginal_survival(d, c(10, 100))
#' td_local_evi(d, 100)
#' @export
transport_dgp <- function(
  x_dist = c("lognormal", "pareto"),
  sdlog = 1,
  alpha = 3,
  xi_cond = 0.2
) {
  x_dist <- match.arg(x_dist)
  if (xi_cond <= 0) {
    stop("xi_cond must be positive.", call. = FALSE)
  }
  if (sdlog <= 0 || alpha <= 0) {
    stop("sdlog and alpha must be positive.", call. = FALSE)
  }
  dgp <- structure(
    list(
      x_dist = x_dist,
      sdlog = sdlog,
      alpha = alpha,
      xi_cond = xi_cond
    ),
    class = "transport_dgp"
  )
  # The quadrature is exact but too slow to sit inside a bisection or a Shiny
  # redraw. log S_Y is smooth and near-linear in log y, so cache it once on a
  # fine log-y grid and interpolate; `exact = TRUE` still goes to the quadrature,
  # and the tests check the two agree.
  evi <- td_asymptotic_evi(dgp)
  ly <- seq(-25, 250 * base::log(10) * evi + 25, length.out = 4000L)
  dgp$cache <- list(
    log_y = ly,
    log_s = td_marginal_survival(dgp, exp(ly), ngrid = 30001L,
                                 log = TRUE, exact = TRUE)
  )
  dgp
}

#' @export
print.transport_dgp <- function(x, ...) {
  cat("<transport_dgp>\n")
  cat(
    "  X ~ ",
    if (x$x_dist == "lognormal") {
      paste0("LogNormal(0, ", format(x$sdlog), ")")
    } else {
      paste0("Pareto(alpha = ", format(x$alpha), ", scale = 1)")
    },
    "\n",
    sep = ""
  )
  cat("  Y | X = x ~ GP(scale = x, shape = ", format(x$xi_cond), ")\n", sep = "")
  cat(
    "  marginal EVI: ",
    format(td_asymptotic_evi(x)),
    " asymptotically, CEVI ",
    format(round(x$xi_cond / td_asymptotic_evi(x), 4)),
    "\n",
    sep = ""
  )
  invisible(x)
}

#' Asymptotic marginal tail index of the DGP
#'
#' Breiman's lemma: `Y = X * W` is regularly varying with index the heavier of
#' the two factors, so the marginal EVI is `max(xi_cond, 1 / alpha)` for a
#' Pareto covariate and `xi_cond` for a lognormal one. The ratio
#' `xi_cond / this` is the copula conditional EVI.
#'
#' @param dgp A `transport_dgp`.
#' @returns A single number.
#' @examples
#' td_asymptotic_evi(transport_dgp("pareto", alpha = 3, xi_cond = 0.2))
#' @export
td_asymptotic_evi <- function(dgp) {
  if (dgp$x_dist == "lognormal") dgp$xi_cond else max(dgp$xi_cond, 1 / dgp$alpha)
}

#' The covariate distribution
#'
#' @inheritParams td_asymptotic_evi
#' @param p Numeric vector of probabilities.
#' @param x Numeric vector of covariate values.
#' @returns A numeric vector.
#' @examples
#' td_x_quantile(transport_dgp(), c(0.1, 0.9))
#' @export
td_x_quantile <- function(dgp, p) {
  if (dgp$x_dist == "lognormal") {
    stats::qlnorm(p, 0, dgp$sdlog)
  } else {
    (1 - p)^(-1 / dgp$alpha)
  }
}

#' @rdname td_x_quantile
#' @export
td_x_cdf <- function(dgp, x) {
  if (dgp$x_dist == "lognormal") {
    stats::plnorm(x, 0, dgp$sdlog)
  } else {
    ifelse(x < 1, 0, 1 - x^(-dgp$alpha))
  }
}

# Log density of L = log X, on the grid the quadrature integrates over.
td_log_density_logx <- function(dgp, l) {
  if (dgp$x_dist == "lognormal") {
    stats::dnorm(l, 0, dgp$sdlog, log = TRUE)
  } else {
    ifelse(l < 0, -Inf, log(dgp$alpha) - dgp$alpha * l)
  }
}

#' Conditional survival of the DGP
#'
#' `P(Y > y | X = x)`, which is exactly generalized Pareto with scale `x`.
#'
#' @inheritParams td_asymptotic_evi
#' @param y,x Numeric vectors, recycled against each other.
#' @returns A numeric vector.
#' @examples
#' td_conditional_survival(transport_dgp(), y = 10, x = 2)
#' @export
td_conditional_survival <- function(dgp, y, x) {
  exp(-(1 / dgp$xi_cond) * log1p(dgp$xi_cond * pmax(y, 0) / x))
}

# Integration grid in l = log x, wide enough for the largest y asked for.
#
# The integrand is exp(-(1/xi) log1p(xi y e^-l) + log f(l)). Below l = log(xi y)
# the first term rises like l / xi and above it flattens, so the mode sits where
# that rise is cancelled by the density's own decay -- at l = sdlog^2 / xi for a
# lognormal, and at l = log(xi y) itself for a Pareto too light to cancel it.
# The grid has to hold both, plus room for the density's decay beyond.
td_grid <- function(dgp, y_max, ngrid) {
  ly <- log(max(dgp$xi_cond * y_max, 1))
  if (dgp$x_dist == "lognormal") {
    pad <- 10 * dgp$sdlog
    lo <- -pad
    hi <- max(pad, dgp$sdlog^2 / dgp$xi_cond + pad, ly + pad)
  } else {
    pad <- 12 / dgp$alpha
    lo <- 0
    hi <- max(pad, ly + pad)
  }
  seq(lo, hi, length.out = ngrid)
}

#' Marginal survival of the DGP
#'
#' `S_Y(y) = E[(1 + xi y / X)^(-1/xi)]`, by quadrature over `log X`. The
#' integration is done in log space with a log-sum-exp, so the result stays
#' accurate down to probabilities far below anything that could be simulated.
#'
#' @inheritParams td_conditional_survival
#' @param ngrid Number of quadrature points.
#' @param log Return `log S_Y(y)` instead of `S_Y(y)`.
#' @param exact Run the quadrature rather than interpolating the curve cached
#'   in the `transport_dgp` object. The cache is accurate to a few parts in a
#'   million and is what makes the app and the bisections usable.
#' @returns A numeric vector.
#' @examples
#' td_marginal_survival(transport_dgp(), c(1, 10, 100))
#' @export
td_marginal_survival <- function(dgp, y, ngrid = 8001L, log = FALSE,
                                 exact = FALSE) {
  y <- pmax(y, 0)
  if (!exact && !is.null(dgp$cache)) {
    out <- stats::approx(
      dgp$cache$log_y,
      dgp$cache$log_s,
      xout = base::log(y),
      rule = 2
    )$y
    # Below the grid the survival is one; above it, extend the last slope.
    below <- base::log(y) < dgp$cache$log_y[1L]
    out[below] <- 0
    n <- length(dgp$cache$log_y)
    above <- base::log(y) > dgp$cache$log_y[n]
    if (any(above)) {
      slope <- diff(dgp$cache$log_s[c(n - 1L, n)]) /
        diff(dgp$cache$log_y[c(n - 1L, n)])
      out[above] <- dgp$cache$log_s[n] +
        slope * (base::log(y[above]) - dgp$cache$log_y[n])
    }
    return(if (log) out else exp(out))
  }
  l <- td_grid(dgp, max(y, 1), ngrid)
  dl <- l[2L] - l[1L]
  lf <- td_log_density_logx(dgp, l)
  keep <- is.finite(lf)
  l <- l[keep]
  lf <- lf[keep]
  eml <- exp(-l)
  # Trapezoid, not rectangle: a Pareto covariate puts a hard edge at log x = 0,
  # and counting that endpoint at full weight costs an O(dl) error there --
  # 3% at the spacings this runs at.
  wt <- rep(dl, length(l))
  wt[c(1L, length(wt))] <- dl / 2
  out <- vapply(
    y,
    function(yy) {
      lg <- -(1 / dgp$xi_cond) * log1p(dgp$xi_cond * yy * eml) + lf
      m <- max(lg)
      if (!is.finite(m)) {
        return(-Inf)
      }
      m + base::log(sum(exp(lg - m) * wt))
    },
    numeric(1L)
  )
  # The quadrature is normalised to the density it integrates, not to one, so a
  # tiny truncation error would show up as S_Y(0) != 1. Divide it out.
  l0 <- {
    m <- max(lf)
    m + base::log(sum(exp(lf - m) * wt))
  }
  out <- out - l0
  if (log) out else exp(out)
}

#' Marginal quantile of the DGP
#'
#' Inverts [td_marginal_survival()] by bisection in `log y`, so it stays exact
#' for exceedance probabilities well beyond the reach of simulation.
#'
#' @inheritParams td_marginal_survival
#' @param s Numeric vector of upper-tail probabilities.
#' @returns A numeric vector.
#' @examples
#' td_marginal_quantile(transport_dgp(), c(1e-2, 1e-6))
#' @export
td_marginal_quantile <- function(dgp, s, ngrid = 8001L, exact = FALSE) {
  if (!exact && !is.null(dgp$cache)) {
    # log_s is strictly decreasing, so inverting it is one interpolation.
    return(exp(stats::approx(
      rev(dgp$cache$log_s),
      rev(dgp$cache$log_y),
      xout = base::log(s),
      rule = 2
    )$y))
  }
  vapply(
    s,
    function(ss) {
      if (!is.finite(ss) || ss <= 0 || ss >= 1) {
        return(NA_real_)
      }
      lo <- -40
      hi <- 40
      for (i in seq_len(200L)) {
        if (td_marginal_survival(dgp, exp(hi), ngrid) > ss) hi <- hi + 20 else break
      }
      for (i in seq_len(200L)) {
        mid <- (lo + hi) / 2
        if (td_marginal_survival(dgp, exp(mid), ngrid) > ss) lo <- mid else hi <- mid
        if (hi - lo < 1e-10) {
          break
        }
      }
      exp((lo + hi) / 2)
    },
    numeric(1L)
  )
}

#' Local tail index of the DGP's marginal
#'
#' `-dlog y / dlog S_Y(y)`, by central difference. This is what a generalized
#' Pareto fitted around `y` would read, and the gap between it and
#' [td_asymptotic_evi()] is the bias built into any canned fit.
#'
#' @inheritParams td_marginal_survival
#' @param eps Half-width of the difference, in `log y`.
#' @returns A numeric vector.
#' @examples
#' td_local_evi(transport_dgp(), c(10, 100, 1e6))
#' @export
td_local_evi <- function(dgp, y, eps = 0.05, ngrid = 8001L) {
  lo <- td_marginal_survival(dgp, y * exp(-eps), ngrid, log = TRUE)
  hi <- td_marginal_survival(dgp, y * exp(eps), ngrid, log = TRUE)
  2 * eps / (lo - hi)
}

#' Simulate from the DGP
#'
#' @inheritParams td_asymptotic_evi
#' @param n Number of observations.
#' @returns A data frame with columns `x`, `y` and `u = F_X(x)`.
#' @examples
#' head(td_simulate(transport_dgp(), 5))
#' @export
td_simulate <- function(dgp, n) {
  u <- stats::runif(n)
  x <- td_x_quantile(dgp, u)
  w <- (stats::runif(n)^(-dgp$xi_cond) - 1) / dgp$xi_cond
  data.frame(x = x, y = x * w, u = u)
}

#' The DGP's own copula, in survival space
#'
#' `td_transport()` takes `P(Y > y | X = x)` and the rank `u = F_X(x)` and
#' returns `P(Y > y)`. It is the exact `h`-inverse of the copula this DGP
#' implies, obtained by undoing the conditional generalized Pareto to recover
#' `y` and then evaluating the derived marginal there -- no copula family is
#' fitted or assumed anywhere. `td_h_surv()` goes the other way.
#'
#' Passing `td_transport` to [marginal_ensemble()] as `family` runs the
#' transport under the true dependence; passing `"gaussian"` or a Clayton family
#' instead runs it under a fitted one. The difference between the two is model
#' risk with everything else held fixed.
#'
#' @inheritParams td_marginal_survival
#' @param s_cond,s_marg Numeric vectors of exceedance probabilities.
#' @param u Numeric vector of `F_X(x)` values.
#' @returns A numeric vector of exceedance probabilities.
#' @examples
#' d <- transport_dgp()
#' td_transport(d, 1e-3, u = 0.9)
#' @export
td_transport <- function(dgp, s_cond, u, ngrid = 8001L) {
  n <- max(length(s_cond), length(u))
  s_cond <- rep_len(s_cond, n)
  u <- rep_len(u, n)
  x <- td_x_quantile(dgp, u)
  # y such that P(Y > y | X = x) = s_cond, in a form that survives tiny s_cond.
  y <- x * expm1(-dgp$xi_cond * base::log(s_cond)) / dgp$xi_cond
  td_marginal_survival(dgp, y, ngrid)
}

#' @rdname td_transport
#' @export
td_h_surv <- function(dgp, s_marg, u, ngrid = 8001L) {
  n <- max(length(s_marg), length(u))
  s_marg <- rep_len(s_marg, n)
  u <- rep_len(u, n)
  y <- td_marginal_quantile(dgp, s_marg, ngrid)
  td_conditional_survival(dgp, y, td_x_quantile(dgp, u))
}

#' Learn conditional tails from a sample, without binning
#'
#' A kernel-weighted conditional distribution at each of `n_anchor` covariate
#' values: every observation gets weight `K((u_i - u_0) / bandwidth)` on the rank
#' scale, and the conditional is the resulting weighted empirical distribution.
#'
#' That is the same object a quantile regression forest returns -- a weighted
#' empirical distribution over the training responses, kernel weights instead of
#' forest weights -- so the tails it produces feed the same machinery the
#' pipeline uses, and there is one distribution per covariate value rather than
#' one per bin.
#'
#' @param u,y Numeric vectors of covariate ranks in `(0, 1)` and responses.
#' @param at Covariate ranks to learn a conditional at. Defaults to
#'   `n_anchor` evenly spaced values.
#' @param n_anchor Number of covariate values, when `at` is not given.
#' @param bandwidth Kernel bandwidth on the rank scale.
#' @param tail_frac Upper fraction of each conditional used for its generalized
#'   Pareto tail.
#' @returns A list with `pieces` (tail pieces in the shape [dl_tail_pieces()]
#'   returns, `NULL` where a conditional is too thin) and `u`.
#' @examples
#' d <- td_simulate(transport_dgp(), 500)
#' length(td_conditionals(d$u, d$y, n_anchor = 10)$pieces)
#' @export
td_conditionals <- function(
  u,
  y,
  at = NULL,
  n_anchor = 40L,
  bandwidth = 0.06,
  tail_frac = 0.35
) {
  if (is.null(at)) {
    at <- seq(0.03, 0.97, length.out = n_anchor)
  }
  o <- order(y)
  ys <- y[o]
  us <- u[o]
  pieces <- lapply(at, function(u0) {
    w <- stats::dnorm((us - u0) / bandwidth)
    tot <- sum(w)
    if (!is.finite(tot) || tot <= 0) {
      return(NULL)
    }
    w <- w / tot
    # Threshold at the tail_frac point of this conditional, read off its own
    # weighted distribution.
    thr <- ys[which(cumsum(w) >= 1 - tail_frac)[1L]]
    above <- ys > thr
    if (sum(above) < 5L) {
      return(NULL)
    }
    wa <- w[above]
    list(
      threshold = thr,
      excess = ys[above] - thr,
      weight = wa,
      tail_prob = sum(wa),
      n_eff = sum(wa)^2 / sum(wa^2)
    )
  })
  list(pieces = pieces, u = at)
}

#' A generalized Pareto fitted to the marginal directly
#'
#' The comparison the transport is against: pool all the responses, ignore the
#' covariate, and fit one generalized Pareto above a high quantile. Under
#' [transport_dgp()] this is exactly the canned fit that cannot work, because
#' the derived marginal is not generalized Pareto anywhere a sample reaches.
#'
#' @param y Numeric vector of responses.
#' @param p Threshold quantile level.
#' @returns A list with `threshold`, `tail_prob`, `scale` and `shape`, and a
#'   `survival(y)` closure.
#' @examples
#' d <- td_simulate(transport_dgp(), 2000)
#' td_canned_gp(d$y)$shape
#' @export
td_canned_gp <- function(y, p = 0.9) {
  y <- y[is.finite(y)]
  thr <- stats::quantile(y, p, names = FALSE)
  e <- y[y > thr] - thr
  fit <- fit_gpd_weighted(e, rep(1, length(e)))
  tail_prob <- length(e) / length(y)
  list(
    threshold = thr,
    tail_prob = tail_prob,
    scale = fit$scale,
    shape = fit$shape,
    survival = function(yy) {
      tail_prob * gpd_survival(yy - thr, fit$scale, fit$shape)
    }
  )
}
