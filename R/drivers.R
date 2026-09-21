#' Joint model of the predictors ("drivers") at peak hours
#'
#' One parametric marginal per predictor ([famish::fit_dst()]), joined by a
#' bivariate copula ([rvinecopulib::bicop()]) when there are two predictors.
#' With one predictor the model is just the marginal. This generalises the
#' earlier rainfall–snowmelt-only joint model to whatever predictors an
#' analysis uses.
#'
#' @param data Training data for one cell.
#' @param predictors Predictor column names (one or two).
#' @param marginal_family Family name passed to [famish::fit_dst()].
#' @param family_set Copula families passed to [rvinecopulib::bicop()].
#' @returns An object of class `"drivers"`: a list with `predictors`,
#'   `marginals` (named list of distributions) and `copula` (a `bicop_dist` or
#'   `NULL`).
#' @export
fit_drivers <- function(data, predictors, marginal_family = "gamma",
                        family_set = "parametric") {
  if (length(predictors) > 2) {
    rlang::abort("fit_drivers() supports one or two predictors for now.")
  }
  data <- tidyr::drop_na(data[predictors])
  marginals <- lapply(predictors, function(p) {
    x <- data[[p]]
    # Zero-inflated predictors: fit the family to a tiny positive offset so
    # families on (0, Inf) still fit; the mass at zero is negligible at peaks.
    famish::fit_dst(marginal_family, pmax(x, 1e-6))
  })
  names(marginals) <- predictors
  copula <- NULL
  if (length(predictors) == 2) {
    u <- cdf_scores_to_copula_u(
      distionary::eval_cdf(marginals[[1]], at = data[[predictors[1]]]),
      distionary::eval_cdf(marginals[[2]], at = data[[predictors[2]]])
    )
    copula <- tryCatch(
      rvinecopulib::bicop(u, family_set = family_set),
      error = function(e) NULL
    )
  }
  structure(
    list(predictors = predictors, marginals = marginals, copula = copula),
    class = c("drivers", "list")
  )
}

#' Evaluate the drivers density
#'
#' @param object A `"drivers"` object.
#' @param newdata Data frame with the predictor columns.
#' @param given For two predictors, optionally the name of the predictor to
#'   condition on: returns the conditional density of the other predictor.
#' @returns Numeric vector of density values.
#' @export
eval_drivers_density <- function(object, newdata, given = NULL) {
  checkmate::assert_class(object, "drivers")
  preds <- object$predictors
  dens <- lapply(preds, function(p) {
    distionary::eval_density(object$marginals[[p]], at = newdata[[p]])
  })
  names(dens) <- preds
  if (length(preds) == 1) {
    return(dens[[1]])
  }
  cop <- 1
  if (!is.null(object$copula)) {
    u <- cdf_scores_to_copula_u(
      distionary::eval_cdf(object$marginals[[1]], at = newdata[[preds[1]]]),
      distionary::eval_cdf(object$marginals[[2]], at = newdata[[preds[2]]])
    )
    bc <- object$copula
    cop <- as.numeric(rvinecopulib::dbicop(u, bc$family, bc$rotation, bc$parameters))
  }
  if (is.null(given)) {
    return(cop * dens[[1]] * dens[[2]])
  }
  checkmate::assert_choice(given, preds)
  cop * dens[[setdiff(preds, given)]]
}

#' Interiorize one margin's CDF scores away from 0 and 1
#'
#' @param u Numeric vector of cdf scores in \eqn{[0, 1]}.
#' @param eps Candidate lower bound for exact zeros.
#' @noRd
interiorize_margin_cdf <- function(u, eps = 1e-6) {
  u <- pmax(pmin(as.numeric(u), 1), 0)
  pos <- u[u > 0]
  lo <- if (length(pos) > 0) min(pos) else eps
  low <- if (eps > lo) lo / 2 else eps
  high <- 1 - low
  u[u == 0] <- low
  u[u == 1] <- high
  u
}

#' Map margin CDF scores into the interior of \eqn{(0,1)} for copula fitting
#'
#' `rvinecopulib::bicop()` expects pseudo-observations strictly inside the unit
#' hypercube. Exact zeros and ones from [distionary::eval_cdf()] are replaced
#' per margin: zeros become `eps` or half the smallest positive score in that
#' margin (whichever is smaller); ones become one minus that lower bound.
#'
#' @param u,v Numeric vectors of margin cdf scores in \eqn{[0, 1]} (same length).
#' @param ... Optional additional margin cdf vectors (same length as `u`).
#' @param eps Lower interior bound for exact zeros when no smaller positive
#'   score exists in a margin.
#' @returns Numeric matrix with one column per margin (`u`, `v`, then `...`).
#' @noRd
cdf_scores_to_copula_u <- function(u, v, ..., eps = 1e-6) {
  rlang::check_dots_empty()
  checkmate::assert_numeric(u)
  checkmate::assert_numeric(v)
  uv <- vctrs::vec_recycle_common(u, v)
  u <- uv[[1]]
  v <- uv[[2]]
  u <- interiorize_margin_cdf(u, eps = eps)
  v <- interiorize_margin_cdf(v, eps = eps)
  cbind(u, v)
}


