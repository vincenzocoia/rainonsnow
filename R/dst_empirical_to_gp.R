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

  init_fit <- famish::fit_dst_gp(excesses)
  par <- distionary::parameters(init_fit)

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
      fit$convergence, ")."
    )
  }

  distionary::dst_gp(scale = exp(fit$par[1]), shape = fit$par[2]) + threshold
}