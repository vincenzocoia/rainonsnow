# Estimators compared at hydrological sample sizes.
#
# Everything here works from a record of n years. `n` is 50 to 100 in hydrology,
# which is the regime that decides whether any of this is useful; larger records
# are included because other fields have them.
#
# The comparison is between three ways of getting a marginal return level:
#
#   1. Fit a distribution to the annual peak directly. Standard practice, and
#      the thing to beat.
#   2. Fit the two drivers separately and recombine. This is the decomposition
#      argument in its original form -- it needs the drivers to be observed.
#   3. Fit the response given a covariate, then get the marginal from that. Two
#      ways to take the last step: integrate the conditional over a FITTED
#      covariate distribution, or transport each conditional through a copula
#      and combine. The first is the law of total probability with the
#      covariate's own tail restored; the second needs no covariate
#      distribution but does need the copula.
#
# The difference between the last two is worth stating plainly, because it is
# the whole reason the mixture in the pipeline is biased. Averaging the
# conditionals over the OBSERVED covariate values truncates the covariate
# distribution at the sample maximum, and everything the marginal tail gets from
# rarer covariate values is lost. Integrating over a fitted covariate
# distribution puts exactly that back, and needs no copula to do it.

#' Generalized extreme value fit and return levels
#'
#' The standard competitor for annual maxima: three parameters from all `n`
#' years, rather than two from the exceedances.
#'
#' @param y Numeric vector of annual maxima.
#' @returns A list with `location`, `scale`, `shape` and `converged`.
#' @examples
#' set.seed(1)
#' fit_gev(stats::rnorm(60, 10, 2))$shape
#' @export
fit_gev <- function(y) {
  y <- y[is.finite(y)]
  nllh <- function(p) {
    mu <- p[1]
    sig <- exp(p[2])
    xi <- p[3]
    z <- (y - mu) / sig
    if (abs(xi) < 1e-8) {
      return(sum(base::log(sig) + z + exp(-z)))
    }
    t <- 1 + xi * z
    if (any(t <= 0)) {
      return(1e10)
    }
    sum(base::log(sig) + (1 + 1 / xi) * base::log(t) + t^(-1 / xi))
  }
  start <- c(stats::median(y), base::log(stats::sd(y) * 0.78), 0.1)
  o <- try(
    stats::optim(start, nllh, control = list(maxit = 2000, reltol = 1e-10)),
    silent = TRUE
  )
  if (inherits(o, "try-error") || o$value >= 1e10) {
    return(list(location = NA_real_, scale = NA_real_, shape = NA_real_,
                converged = FALSE))
  }
  list(
    location = o$par[1],
    scale = exp(o$par[2]),
    shape = o$par[3],
    converged = o$convergence == 0
  )
}

#' @rdname fit_gev
#' @param fit A `fit_gev()` result.
#' @param p Numeric vector of annual exceedance probabilities (`1 / T`).
#' @export
gev_return_level <- function(fit, p) {
  if (!is.finite(fit$shape)) {
    return(rep(NA_real_, length(p)))
  }
  yp <- -base::log(1 - p)
  if (abs(fit$shape) < 1e-8) {
    return(fit$location - fit$scale * base::log(yp))
  }
  fit$location + fit$scale / fit$shape * (yp^(-fit$shape) - 1)
}

#' Generalized Pareto fit to the whole of a positive sample
#'
#' Used for the rainfall marginal, whose tail is what the covariate
#' distribution has to supply when the conditional is integrated over it.
#'
#' @param x Numeric vector of positive values.
#' @returns A list with `scale` and `shape`.
#' @examples
#' set.seed(1)
#' fit_gp_ml(stats::rexp(200))$shape
#' @export
fit_gp_ml <- function(x) {
  x <- x[is.finite(x) & x > 0]
  nllh <- function(p) {
    d <- gpd_density(x, exp(p[1]), p[2])
    if (any(!is.finite(d)) || any(d <= 0)) {
      return(1e10)
    }
    -sum(base::log(d))
  }
  o <- try(
    stats::optim(c(base::log(stats::median(x)), 0.1), nllh,
                 control = list(maxit = 2000, reltol = 1e-10)),
    silent = TRUE
  )
  if (inherits(o, "try-error") || o$value >= 1e10) {
    return(list(scale = NA_real_, shape = NA_real_))
  }
  list(scale = exp(o$par[1]), shape = o$par[2])
}

#' Fit the two-driver structure to peaks and a covariate
#'
#' Assumes what [two_process_dgp()] generates: an unobserved driver `A` that
#' does not depend on the covariate, a driver `B` whose scale is proportional to
#' it, and a peak equal to their maximum. Four parameters, fitted by maximum
#' likelihood from `(x, y)` pairs alone -- neither driver is observed.
#'
#' The two drivers are exchangeable without a constraint, so the fit imposes
#' `shape_a <= shape_b`: `A` is the lighter one. In the intended application it
#' is also the bounded one, since snowmelt is limited by the snowpack.
#'
#' @param x,y Numeric vectors: the covariate and the annual peak.
#' @returns A list with `scale_a`, `shape_a`, `b_coef`, `shape_b`, `nllh` and
#'   `converged`.
#' @examples
#' set.seed(1)
#' d <- tp_simulate(two_process_dgp(), 400)
#' fit_two_process(d$x, d$y)$shape_b
#' @export
fit_two_process <- function(x, y) {
  ok <- is.finite(x) & is.finite(y) & x > 0 & y > 0
  x <- x[ok]
  y <- y[ok]
  unpack <- function(p) {
    list(
      scale_a = exp(p[1]),
      shape_a = p[2],
      b_coef = exp(p[3]),
      # shape_b >= shape_a, so the drivers cannot swap labels.
      shape_b = p[2] + base::log1p(exp(p[4]))
    )
  }
  nllh <- function(p) {
    q <- unpack(p)
    sa <- gpd_survival(y, q$scale_a, q$shape_a)
    sb <- gpd_survival(y, q$b_coef * x, q$shape_b)
    da <- gpd_density(y, q$scale_a, q$shape_a)
    db <- gpd_density(y, q$b_coef * x, q$shape_b)
    dens <- da * (1 - sb) + (1 - sa) * db
    if (any(!is.finite(dens)) || any(dens <= 0)) {
      return(1e10)
    }
    -sum(base::log(dens))
  }
  best <- NULL
  for (s0 in list(c(0, -0.2, -0.7, 0), c(0.5, -0.1, -1.5, -1), c(-0.5, 0.05, 0, 1))) {
    start <- c(base::log(stats::median(y)) + s0[1], s0[2], s0[3], s0[4])
    o <- try(
      stats::optim(start, nllh, control = list(maxit = 4000, reltol = 1e-11)),
      silent = TRUE
    )
    if (!inherits(o, "try-error") && o$value < 1e10 &&
        (is.null(best) || o$value < best$value)) {
      best <- o
    }
  }
  if (is.null(best)) {
    return(list(scale_a = NA_real_, shape_a = NA_real_, b_coef = NA_real_,
                shape_b = NA_real_, nllh = Inf, converged = FALSE))
  }
  c(unpack(best$par), list(nllh = best$value, converged = best$convergence == 0))
}

#' Conditional survival implied by a fitted two-driver model
#'
#' @param fit A [fit_two_process()] result.
#' @param y,x Numeric vectors, recycled against each other.
#' @returns A numeric vector.
#' @examples
#' set.seed(1)
#' d <- tp_simulate(two_process_dgp(), 300)
#' tpf_conditional_survival(fit_two_process(d$x, d$y), y = 5, x = 1)
#' @export
tpf_conditional_survival <- function(fit, y, x) {
  sa <- gpd_survival(y, fit$scale_a, fit$shape_a)
  sb <- gpd_survival(y, fit$b_coef * x, fit$shape_b)
  sa + sb - sa * sb
}

#' Marginal survival, two ways of removing the covariate
#'
#' `tpf_marginal_survival()` integrates the fitted conditional over a FITTED
#' covariate distribution, so the covariate's own tail contributes. This is the
#' step the pipeline's mixture is missing.
#'
#' `tpf_mixture_survival()` averages over the OBSERVED covariate values, which
#' is what the pipeline does now. It truncates the covariate distribution at the
#' sample maximum, and everything the marginal tail owes to rarer covariate
#' values is lost -- a bias that does not shrink with more data, because more
#' data does not reach beyond its own maximum by any margin that matters.
#'
#' @inheritParams tpf_conditional_survival
#' @param fit_x A [fit_gp_ml()] result for the covariate.
#' @param ngrid Quadrature points.
#' @param x_obs Observed covariate values.
#' @returns A numeric vector.
#' @examples
#' set.seed(1)
#' d <- tp_simulate(two_process_dgp(), 300)
#' f <- fit_two_process(d$x, d$y)
#' tpf_marginal_survival(f, fit_gp_ml(d$x), y = c(5, 10))
#' @export
tpf_marginal_survival <- function(fit, fit_x, y, ngrid = 4001L) {
  if (!is.finite(fit$shape_b) || !is.finite(fit_x$shape)) {
    return(rep(NA_real_, length(y)))
  }
  l <- seq(-14, 50, length.out = ngrid)
  dl <- l[2L] - l[1L]
  lf <- base::log(gpd_density(exp(l), fit_x$scale, fit_x$shape)) + l
  ok <- is.finite(lf)
  l <- l[ok]
  w <- rep(dl, length(l))
  w[c(1L, length(w))] <- dl / 2
  fx <- exp(lf[ok]) * w
  total <- sum(fx)
  if (!is.finite(total) || total <= 0) {
    return(rep(NA_real_, length(y)))
  }
  x <- exp(l)
  vapply(
    y,
    function(yy) sum(fx * tpf_conditional_survival(fit, yy, x)) / total,
    numeric(1L)
  )
}

#' @rdname tpf_marginal_survival
#' @export
tpf_mixture_survival <- function(fit, x_obs, y) {
  vapply(
    y,
    function(yy) mean(tpf_conditional_survival(fit, yy, x_obs)),
    numeric(1L)
  )
}

#' Invert a survival function for return levels
#'
#' @param s_fun A function of `y` returning an exceedance probability.
#' @param p Numeric vector of exceedance probabilities.
#' @param upper Largest value to search.
#' @returns A numeric vector, `NA` where `p` is unreachable.
#' @examples
#' invert_survival(function(y) exp(-y), c(0.05, 0.01))
#' @export
invert_survival <- function(s_fun, p, upper = 1e9) {
  vapply(
    p,
    function(pp) {
      lo <- 0
      hi <- upper
      if (!is.finite(s_fun(hi)) || s_fun(hi) > pp) {
        return(NA_real_)
      }
      for (i in seq_len(200L)) {
        m <- (lo + hi) / 2
        s <- s_fun(m)
        if (!is.finite(s)) {
          return(NA_real_)
        }
        if (s > pp) lo <- m else hi <- m
        if (hi - lo < 1e-8 * max(hi, 1)) {
          break
        }
      }
      (lo + hi) / 2
    },
    numeric(1L)
  )
}

#' Fit the rain-driven component from the rainfall that caused it
#'
#' The marginal distribution of a driver whose scale is set by a random
#' covariate is itself a scale mixture, so fitting a generalized Pareto to it
#' is misspecified even though the conditional one is exact. This fits the
#' conditional: scale proportional to the covariate, one shape.
#'
#' @param x,b Numeric vectors: the covariate and the driver it scales.
#' @returns A list with `b_coef` and `shape`.
#' @examples
#' set.seed(1)
#' d <- tp_simulate(two_process_dgp(), 400)
#' fit_scaled_gp(d$x, d$b)$b_coef
#' @export
fit_scaled_gp <- function(x, b) {
  ok <- is.finite(x) & is.finite(b) & x > 0 & b > 0
  x <- x[ok]
  b <- b[ok]
  nllh <- function(p) {
    d <- gpd_density(b, exp(p[1]) * x, p[2])
    if (any(!is.finite(d)) || any(d <= 0)) {
      return(1e10)
    }
    -sum(base::log(d))
  }
  o <- try(
    stats::optim(c(base::log(stats::median(b / x)), 0.1), nllh,
                 control = list(maxit = 2000, reltol = 1e-10)),
    silent = TRUE
  )
  if (inherits(o, "try-error") || o$value >= 1e10) {
    return(list(b_coef = NA_real_, shape = NA_real_))
  }
  list(b_coef = exp(o$par[1]), shape = o$par[2])
}
