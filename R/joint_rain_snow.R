#' Map margin CDF values into the interior of \eqn{[0,1]} for copula fitting
#'
#' `rvinecopulib::bicop()` expects pseudo-observations strictly inside the unit
#' hypercube. Semi-parametric margins often yield exact zeros and ones; scaling
#' by \eqn{n/(n+1)} matches the usual rank-based pseudo-observation adjustment.
#'
#' @param u_mat Numeric matrix with two columns (rain, snowmelt margin cdf
#'   scores).
#' @param n Sample size used for scaling.
#' @returns A numeric matrix of the same shape as `u_mat`.
#' @noRd
cdf_scores_to_copula_u <- function(u_mat, n) {
  u_mat <- pmax(pmin(as.matrix(u_mat), 1), 0)
  u_mat * (n / (n + 1L))
}


#' Joint model for hourly rainfall and snowmelt at one grid cell
#'
#' Fits two **famish** marginals (each a `distionary` distribution object) and links
#' them with a bivariate copula from **rvinecopulib** (`bicop()`). Because the
#' probaverse stack does not yet represent multivariate distributions as a
#' single `dst`, the result is a bespoke list with S3 class `joint_rain_snow`.
#'
#' Pseudo-observations are margin cdf scores from `distionary::eval_cdf()`,
#' scaled to fall strictly inside \eqn{(0,1)^2} so `bicop()` accepts them.
#'
#' @param rainfall_hourly,snowmelt_hourly Paired numeric vectors (same length),
#'   typically mm/h as in `scripts/2-tablify_spatial_eo.r`.
#' @param cell_id,x,y Optional coordinates / id stored in `meta` for book-keeping.
#' @param marginal_rainfall,marginal_snowmelt Distribution **family** names passed to
#'   [famish::fit_dst()] (e.g. `"empirical"`, `"gamma"`, `"weibull"`).
#' @param bicop_family_set,bicop_controls Passed to [rvinecopulib::bicop()].
#' @param min_obs Minimum paired complete observations required to fit the
#'   copula; otherwise `bicop` is `NULL`.
#' @returns An object of class `joint_rain_snow`: a list with components
#'   `marginal_rainfall`, `marginal_snowmelt`, `bicop` (may be `NULL`),
#'   `meta` (includes `n_obs`, families, and optional cell fields).
#'
#' @seealso [simulate_joint_rain_snow]
#' @export
fit_joint_rain_snow_cell <- function(
    rainfall_hourly,
    snowmelt_hourly,
    cell_id = NA_integer_,
    x = NA_real_,
    y = NA_real_,
    marginal_rainfall = "empirical",
    marginal_snowmelt = "empirical",
    bicop_family_set = "parametric",
    bicop_controls = list(),
    min_obs = 40L) {
  checkmate::assert_numeric(rainfall_hourly)
  checkmate::assert_numeric(snowmelt_hourly)
  checkmate::assert_true(length(rainfall_hourly) == length(snowmelt_hourly))

  checkmate::assert_string(marginal_rainfall)
  checkmate::assert_string(marginal_snowmelt)

  ok <- stats::complete.cases(rainfall_hourly, snowmelt_hourly)
  rain <- rainfall_hourly[ok]
  sm <- snowmelt_hourly[ok]
  n <- length(rain)

  dst_rain <- famish::fit_dst(marginal_rainfall, rain)
  dst_sm <- famish::fit_dst(marginal_snowmelt, sm)

  meta <- list(
    n_obs = n,
    marginal_rainfall = marginal_rainfall,
    marginal_snowmelt = marginal_snowmelt,
    cell_id = cell_id,
    x = x,
    y = y
  )

  if (n < min_obs) {
    warning(
      "Only ", n, " paired observations (minimum ", min_obs, "); skipping copula."
    )
    return(structure(
      list(
        marginal_rainfall = dst_rain,
        marginal_snowmelt = dst_sm,
        bicop = NULL,
        meta = meta
      ),
      class = c("joint_rain_snow", "list")
    ))
  }

  u1 <- distionary::eval_cdf(dst_rain, at = rain)
  u2 <- distionary::eval_cdf(dst_sm, at = sm)
  u <- cdf_scores_to_copula_u(cbind(u1, u2), n)

  bc <- tryCatch(
    do.call(
      rvinecopulib::bicop,
      c(
        list(data = u, family_set = bicop_family_set),
        bicop_controls
      )
    ),
    error = function(e) {
      warning("Could not fit bivariate copula: ", conditionMessage(e))
      NULL
    }
  )

  structure(
    list(
      marginal_rainfall = dst_rain,
      marginal_snowmelt = dst_sm,
      bicop = bc,
      meta = meta
    ),
    class = c("joint_rain_snow", "list")
  )
}


#' Fit joint rainfall–snowmelt models for every spatial group in a table
#'
#' Groups `data` (usually `cell_id` or `(x, y)`), fits one `joint_rain_snow`
#' model per group using `fit_joint_rain_snow_cell()`, and returns a row per
#' group with a list column `joint`.
#'
#' @param data A data frame containing paired hourly columns.
#' @param group_cols Character vector of grouping column names (`cell_id`,
#'   `x`, `y`, etc.).
#' @param rainfall_col,snowmelt_col Names of the rainfall and snowmelt columns.
#' @inheritParams fit_joint_rain_snow_cell
#'
#' @param progress If `TRUE`, shows a terminal progress bar over grid cells (via
#'   [utils::txtProgressBar()]).
#' @param verbose If `TRUE`, prints a short [message()] after each cell is
#'   fitted (can be very chatty with many cells).
#'
#' @returns A grouped-by-row data frame with the same grouping keys plus column
#'   `joint` (list of `joint_rain_snow` objects).
#' @export
fit_joint_rain_snow_cells <- function(
    data,
    group_cols,
    rainfall_col = "rainfall_hourly",
    snowmelt_col = "snowmelt_hourly",
    marginal_rainfall = "empirical",
    marginal_snowmelt = "empirical",
    bicop_family_set = "parametric",
    bicop_controls = list(),
    min_obs = 40L,
    progress = FALSE,
    verbose = FALSE) {
  checkmate::assert_data_frame(data)
  checkmate::assert_character(group_cols, min.len = 1L)
  checkmate::assert_subset(group_cols, names(data))
  checkmate::assert_subset(c(rainfall_col, snowmelt_col), names(data))
  checkmate::assert_string(marginal_rainfall)
  checkmate::assert_string(marginal_snowmelt)
  checkmate::assert_flag(progress)
  checkmate::assert_flag(verbose)

  grp <- dplyr::group_by(data, dplyr::across(dplyr::all_of(group_cols)))
  keys <- dplyr::group_keys(grp)
  chunks <- dplyr::group_split(grp)
  n_grp <- nrow(keys)

  if (n_grp == 0L) {
    out <- dplyr::mutate(keys, joint = vector("list", 0L))
    return(dplyr::group_by(
      out,
      dplyr::across(dplyr::all_of(group_cols))
    ))
  }

  if (progress && n_grp > 0L) {
    pb <- utils::txtProgressBar(min = 0L, max = n_grp, style = 3L)
    on.exit(close(pb), add = TRUE)
  }

  rows <- vector("list", n_grp)
  for (i in seq_len(n_grp)) {
    g <- chunks[[i]]
    key <- keys[i, , drop = FALSE]

    cid <- if ("cell_id" %in% names(key)) {
      key[["cell_id"]][[1]]
    } else if ("cell_id" %in% names(g)) {
      g[["cell_id"]][[1]]
    } else {
      NA_integer_
    }
    gx <- if ("x" %in% names(key)) {
      key[["x"]][[1]]
    } else if ("x" %in% names(g)) {
      g[["x"]][[1]]
    } else {
      NA_real_
    }
    gy <- if ("y" %in% names(key)) {
      key[["y"]][[1]]
    } else if ("y" %in% names(g)) {
      g[["y"]][[1]]
    } else {
      NA_real_
    }

    j <- fit_joint_rain_snow_cell(
      g[[rainfall_col]],
      g[[snowmelt_col]],
      cell_id = cid,
      x = gx,
      y = gy,
      marginal_rainfall = marginal_rainfall,
      marginal_snowmelt = marginal_snowmelt,
      bicop_family_set = bicop_family_set,
      bicop_controls = bicop_controls,
      min_obs = min_obs
    )

    rows[[i]] <- data.frame(joint = I(list(j)))

    if (verbose) {
      message(
        sprintf(
          "joint_rain_snow: cell %d / %d (cell_id=%s)",
          i,
          n_grp,
          if (is.na(cid)) "NA" else as.character(cid)
        )
      )
    }
    if (progress && n_grp > 0L) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  if (progress && n_grp > 0L) {
    cat("\n")
  }

  out <- dplyr::bind_cols(keys, dplyr::bind_rows(rows))
  dplyr::group_by(out, dplyr::across(dplyr::all_of(group_cols)))
}


#' Simulate hourly rainfall and snowmelt from a fitted joint model
#'
#' Draws uniform pairs from `joint$bicop` with [rvinecopulib::rbicop()] and maps
#' them through the inverse cdf of each **famish** marginal.
#'
#' @param joint An object from `fit_joint_rain_snow_cell()` or from the `joint`
#'   column after `fit_joint_rain_snow_cells()`.
#' @param n Number of simulated pairs.
#' @returns A data frame with columns `rainfall_hourly` and `snowmelt_hourly`.
#' @export
simulate_joint_rain_snow <- function(joint, n) {
  checkmate::assert_class(joint, "joint_rain_snow")
  checkmate::assert_count(n, positive = TRUE)
  if (is.null(joint$bicop)) {
    stop("No copula was fitted (`bicop` is NULL); cannot simulate.", call. = FALSE)
  }

  bc <- joint$bicop
  uv <- rvinecopulib::rbicop(n, bc$family, bc$rotation, bc$parameters)
  rain <- distionary::eval_quantile(joint$marginal_rainfall, at = uv[, 1])
  sm <- distionary::eval_quantile(joint$marginal_snowmelt, at = uv[, 2])
  data.frame(
    rainfall_hourly = rain,
    snowmelt_hourly = sm
  )
}


#' Print a `joint_rain_snow` object
#'
#' @param x A `joint_rain_snow` object from [fit_joint_rain_snow_cell()].
#' @param ... Unused.
#'
#' @returns `x`, invisibly.
#' @export
print.joint_rain_snow <- function(x, ...) {
  cat("<joint_rain_snow>\n")
  cat("  n_obs: ", x$meta$n_obs, "\n", sep = "")
  cat("  marginal rainfall: ", x$meta$marginal_rainfall, "\n", sep = "")
  cat("  marginal snowmelt: ", x$meta$marginal_snowmelt, "\n", sep = "")
  if (!is.null(x$bicop)) {
    cat("  bicop family: ", x$bicop$family, "\n", sep = "")
  } else {
    cat("  bicop: not fitted\n")
  }
  invisible(x)
}
