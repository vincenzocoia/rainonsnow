#' Local linear quantile regression distributional learning model
#'
#' Fits a conditional distribution by kernel-weighted local linear quantile
#' regression, evaluated across a grid of quantile levels. This is the
#' "LOESS, but with quantile regression instead of least squares" idea: for a
#' given `x0`, weight the observations by their distance from `x0`, fit a
#' straight line through the weighted data by minimising the pinball loss at
#' level `p`, and read the fitted line off at `x0`. Repeating over a grid of
#' `p` traces out the whole conditional quantile function, which is then
#' returned as a step distribution.
#'
#' The method is that of Yu and Jones (1998), who localise the characterisation
#' of a regression quantile as the minimiser of `E{rho_p(Y - a) | X = x}`. It
#' is used here as a distributional learning model: unlike a random forest, the
#' conditional quantile function is a smooth function of `x` within the span,
#' so neighbouring `x` give neighbouring distributions.
#'
#' Fitted quantile curves at different levels can cross. Crossings are removed
#' by monotone rearrangement — sorting the fitted quantiles at each `x0` — which
#' cannot increase the estimation error (Chernozhukov, Fernandez-Val and
#' Galichon, 2010).
#'
#' Only one predictor is supported.
#'
#' @param data A data frame.
#' @param yname The name of the response variable.
#' @param xnames The name of the (single) predictor variable.
#' @param levels Quantile levels at which to fit. The returned distribution is
#'   a step function with these levels as its jump points.
#' @param span Fraction of the data falling inside the kernel window, as in
#'   `stats::loess()`. Ignored if `bandwidth` is given.
#' @param bandwidth Kernel half-width in the units of `x`. If `NULL` (the
#'   default) a nearest-neighbour bandwidth is used at each `x0`, giving the
#'   locally adaptive window that `span` describes.
#' @param kernel Kernel used to weight observations by distance from `x0`.
#' @param degree Local polynomial degree: 1 for local linear (the default), or
#'   0 for a kernel-weighted moving window with no slope term.
#' @param na_action Either "drop" to drop rows with missing values, or "null"
#'   to return a "dl_null" object.
#' @param min_obs The minimum number of observations required to fit.
#' @returns A distributional learning model with subclass "dl_llqr".
#' @references
#' Yu, K. and Jones, M. C. (1998). Local linear quantile regression.
#' *Journal of the American Statistical Association*, 93(441), 228-237.
#'
#' Chernozhukov, V., Fernandez-Val, I. and Galichon, A. (2010). Quantile and
#' probability curves without crossing. *Econometrica*, 78(3), 1093-1125.
#' @examples
#' set.seed(1)
#' df <- data.frame(x = runif(200, 0, 10))
#' df$y <- df$x + rexp(200, rate = 1 / (1 + df$x / 5))
#' fit <- dl_llqr(df, yname = "y", xnames = "x")
#' predict(fit, newdata = data.frame(x = c(2, 8)))
#' @export
dl_llqr <- function(data, yname, xnames,
                    levels = seq(0.01, 0.99, by = 0.01),
                    span = 0.4,
                    bandwidth = NULL,
                    kernel = c("tricube", "epanechnikov", "gaussian"),
                    degree = 1L,
                    na_action = c("drop", "null"),
                    min_obs = 20) {
  kernel <- rlang::arg_match(kernel)
  na_action <- rlang::arg_match(na_action)
  checkmate::assert_number(span, lower = 0.01, upper = 1)
  checkmate::assert_numeric(levels, lower = 0, upper = 1, any.missing = FALSE)
  checkmate::assert_int(degree, lower = 0, upper = 1)
  if (length(xnames) != 1L) {
    rlang::abort("`dl_llqr()` supports exactly one predictor.")
  }
  data <- data[append(xnames, yname)]
  data2 <- tidyr::drop_na(data)
  if (na_action == "null" && nrow(data2) < nrow(data)) {
    return(dl_null())
  }
  data <- data2
  if (nrow(data) < min_obs) {
    return(dl_null(data))
  }
  res <- list(
    yname = yname,
    xnames = xnames,
    levels = sort(unique(levels)),
    span = span,
    bandwidth = bandwidth,
    kernel = kernel,
    degree = as.integer(degree),
    training = data
  )
  new_dstlrn(res, subclass = "dl_llqr")
}

#' Kernel weights for local fitting
#'
#' @param x Numeric vector of predictor values.
#' @param x0 Scalar point at which to centre the kernel.
#' @param h Kernel half-width. Observations further than `h` from `x0` get zero
#'   weight for the compactly supported kernels.
#' @param kernel Kernel name.
#' @returns A numeric vector of weights the same length as `x`.
#' @examples
#' llqr_weights(c(1, 2, 3), x0 = 2, h = 1.5, kernel = "tricube")
#' @export
llqr_weights <- function(x, x0, h, kernel = "tricube") {
  u <- abs(x - x0) / h
  switch(kernel,
    tricube = ifelse(u < 1, (1 - u^3)^3, 0),
    epanechnikov = ifelse(u < 1, 0.75 * (1 - u^2), 0),
    gaussian = stats::dnorm(u),
    rlang::abort("Unknown kernel.")
  )
}

# Nearest-neighbour bandwidth: the distance to the ceiling(span * n)-th closest
# point, as stats::loess() uses. Guarantees a workable window everywhere,
# including in the sparse right tail of x where a fixed bandwidth would empty.
llqr_nn_bandwidth <- function(x, x0, span) {
  d <- sort(abs(x - x0))
  k <- max(3L, ceiling(span * length(x)))
  d[min(k, length(d))] * 1.0001
}

# Weighted local polynomial quantile regression at a single x0, for a whole
# vector of levels, by the MM algorithm of Hunter and Lange (2000): the pinball
# loss is majorised by a quadratic, so each iteration is a weighted least
# squares solve. Two (or one) parameters, so each solve is closed form.
#
# Minimise  sum_i w_i rho_p(y_i - a - b (x_i - x0))  over (a, b); return a.
llqr_solve_one <- function(xc, y, w, levels, degree,
                           tol = 1e-7, maxit = 60L) {
  keep <- w > 0
  xc <- xc[keep]; y <- y[keep]; w <- w[keep]
  n <- length(y)
  if (n < 3L) return(rep(NA_real_, length(levels)))
  X <- if (degree == 1L) cbind(1, xc) else cbind(rep(1, n))
  Xw <- X * w
  cw <- colSums(Xw)                       # X'w, reused at every level
  eps0 <- 1e-4 * stats::sd(y)
  if (!is.finite(eps0) || eps0 <= 0) eps0 <- 1e-6
  # Start from the weighted least squares fit; refine per level.
  beta_ls <- tryCatch(
    solve(crossprod(X, Xw), crossprod(Xw, y)),
    error = function(e) NULL
  )
  if (is.null(beta_ls)) return(rep(NA_real_, length(levels)))
  out <- numeric(length(levels))
  for (j in seq_along(levels)) {
    p <- levels[j]
    beta <- beta_ls
    eps <- eps0
    for (it in seq_len(maxit)) {
      r <- as.numeric(y - X %*% beta)
      v <- w / (eps + abs(r))
      Xv <- X * v
      A <- crossprod(X, Xv)
      b <- crossprod(Xv, y) - (1 - 2 * p) * cw
      beta_new <- tryCatch(solve(A, b), error = function(e) NULL)
      if (is.null(beta_new)) break
      delta <- max(abs(beta_new - beta))
      beta <- beta_new
      eps <- max(eps * 0.7, 1e-9)
      if (delta < tol) break
    }
    out[j] <- beta[1]                     # value of the local fit at x0
  }
  out
}

#' Conditional quantiles from a local linear quantile regression fit
#'
#' Returns the fitted conditional quantile function evaluated on a grid,
#' before it is turned into a distribution. Useful for plotting the fitted
#' quantile curves themselves.
#'
#' @param object A `"dl_llqr"` object.
#' @param x0 Numeric vector of predictor values at which to evaluate.
#' @param levels Quantile levels; defaults to those stored in `object`.
#' @param rearrange Sort the fitted quantiles at each `x0` to remove crossings.
#' @returns A matrix with one row per element of `x0` and one column per level.
#' @examples
#' set.seed(1)
#' df <- data.frame(x = runif(100, 0, 10))
#' df$y <- df$x + rnorm(100)
#' fit <- dl_llqr(df, yname = "y", xnames = "x", levels = c(0.25, 0.5, 0.75))
#' llqr_quantiles(fit, x0 = c(3, 7))
#' @export
llqr_quantiles <- function(object, x0, levels = NULL, rearrange = TRUE) {
  checkmate::assert_class(object, "dl_llqr")
  levels <- levels %||% object[["levels"]]
  x <- object[["training"]][[object[["xnames"]]]]
  y <- object[["training"]][[object[["yname"]]]]
  h_fixed <- object[["bandwidth"]]
  res <- matrix(NA_real_, nrow = length(x0), ncol = length(levels))
  for (i in seq_along(x0)) {
    h <- h_fixed %||% llqr_nn_bandwidth(x, x0[i], object[["span"]])
    w <- llqr_weights(x, x0[i], h, object[["kernel"]])
    q <- llqr_solve_one(x - x0[i], y, w, levels, object[["degree"]])
    if (rearrange && !anyNA(q)) q <- sort(q)
    res[i, ] <- q
  }
  dimnames(res) <- list(NULL, format(levels))
  res
}

#' @describeIn predict.dstlrn Predict from a local linear quantile regression
#'   distributional learning model. Each prediction is a step distribution
#'   whose jumps sit at the fitted conditional quantiles.
#' @export
predict.dl_llqr <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    newdata <- object[["training"]]
  }
  xname <- object[["xnames"]]
  newdata <- newdata[xname]
  n <- nrow(newdata)
  res <- rep(list(distionary::dst_null()), n)
  lgl_na <- df_rows_have_missing(newdata)
  if (all(lgl_na)) {
    return(res)
  }
  x0 <- newdata[[xname]][!lgl_na]
  qmat <- llqr_quantiles(object, x0)
  levels <- object[["levels"]]
  # Turn the level grid into probability mass: a quantile at level p carries
  # the mass between it and the previous level, with the first knot taking all
  # the mass below it. The result is the step survival function that the
  # tail-grafting machinery expects.
  steps <- diff(c(0, levels))
  steps[length(steps)] <- steps[length(steps)] + (1 - levels[length(levels)])
  dsts <- lapply(seq_len(nrow(qmat)), function(i) {
    q <- qmat[i, ]
    if (anyNA(q)) return(distionary::dst_null())
    distionary::dst_empirical(q, weights = steps)
  })
  res[!lgl_na] <- dsts
  res
}
