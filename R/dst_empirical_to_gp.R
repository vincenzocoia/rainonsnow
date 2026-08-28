#' Empirical distribution to generalized Pareto (excess model)
#'
#' Fits a generalized Pareto distribution to the upper tail of an empirical
#' distribution produced by [distionary::dst_empirical()], using weighted
#' maximum likelihood on the discrete exceedance probabilities. The returned
#' [distionary::dst_gp()] object describes **excesses** `Y = X - u` conditional
#' on `X > u`, with support starting at `0`, matching `dst_gp`'s usual peak-over
#' threshold parameterization.
#'
#' @param dst An empirical distribution: the object returned by
#'   [distionary::dst_empirical()].
#' @param threshold Optional fixed threshold `u`. Excesses are `outcomes - u`
#'   for outcomes strictly greater than `u`. If `NULL`, `u` is chosen
#'   automatically (see `min_tail_prob`).
#' @param min_tail_prob When `threshold` is `NULL`, the largest value `u` on
#'   the support of the empirical distribution such that `P(X > u) >=
#'   min_tail_prob`. Smaller values keep a larger tail sample; larger values
#'   force a higher threshold.
#'
#' @returns A [distionary::dst_gp()] distribution for the fitted scale and shape
#'   (excesses above `u`). Attributes:
#'   * `threshold` — numeric, the threshold `u` used.
#'   * `tail_prob` — `sum(probs[outcomes > u])`, the empirical mass in the tail.
#'   * `n_outcomes` — number of distinct outcome values strictly above `u`.
#'
#' @examples
#' set.seed(1)
#' y <- c(stats::rnorm(200), stats::rnorm(100, mean = 6, sd = 1.2))
#' emp <- distionary::dst_empirical(y)
#' gp <- dst_empirical_to_gp(emp, threshold = stats::quantile(y, 0.85))
#' attr(gp, "threshold")
#' gp
#' @export
convert_emp_to_gp <- function(dst) {
  checkmate::assert_class(dst, "dst")

  pars <- distionary::parameters(dst)
  outcomes <- pars$outcomes
  probs <- pars$probs

  threshold <- range(dst)[1]
  excesses <- outcomes - threshold

  w <- probs

  init_fit <- famish::fit_dst_gp(excesses, method = "mle")

  # Sometimes the excesses go above the max value of the GPD
  # by some machine tolerance.
  rng <- range(init_fit)
  mx <- rng[2]
  tol <- 1e-9
  close <- abs(mx - excesses) < tol
  excesses[close] <- mx - tol

  par <- distionary::parameters(init_fit)
  par <- par[c("scale", "shape")]
  par$scale <- log(par$scale)
  par <- unname(unlist(par))

  negll <- function(par) {
    sigma <- exp(par[1])
    xi <- par[2]
    d <- distionary::dgp(excesses, scale = sigma, shape = xi)
    -sum(w * log(d))
  }

  fit <- stats::optim(par, negll)

  if (fit$convergence != 0L) {
    warning(
      "`optim` reported non-zero convergence code when fitting the GPD (code ",
      fit$convergence,
      ")."
    )
  }

  distionary::dst_gp(scale = exp(fit$par[1]), shape = fit$par[2]) + threshold
}

#' Fit and Graft a GPD to the upper tail of an empirical distribution.
#'
#' @param dst An empirical distribution: the object returned by
#'   `distionary::dst_empirical()`.
#' @param adaptive_threshold Number between 0 and 1. Takes the smallest
#'   value between the quantile of `dst` (with the built-in weights)
#'   and the quantile as if all weights were equal. (This prevents
#'   situations where the empirical tail is dominated by a few outcomes.)
#'
#' @returns A distionary distribution with the GPD grafted back on.
#' @note Requires distplyr > 0.2.0; may have to use the development version.
#' @export
fit_and_graft_gp <- function(dst, adaptive_threshold) {
  checkmate::assert_class(dst, "dst")
  checkmate::assert_number(adaptive_threshold, lower = 0, upper = 1)
  checkmate::assert_true(distionary::pretty_name(dst) == "Finite")
  trim_point1 <- distionary::eval_quantile(dst, at = adaptive_threshold)
  parms <- distionary::parameters(dst)
  outs <- parms$outcomes
  trim_point2 <- quantile(outs, adaptive_threshold)
  trim_point <- min(trim_point1, trim_point2)
  upper_empirical <- distplyr::trim_left(
    dst,
    of = trim_point,
    include = FALSE
  )
  gp <- convert_emp_to_gp(upper_empirical)
  distplyr::graft_right(dst, gp, threshold = trim_point, include = FALSE)
}
