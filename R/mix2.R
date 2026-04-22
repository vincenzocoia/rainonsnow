# Scalar bound from range(dst): j = 1 (lower) or 2 (upper).
# range(dst) can be a numeric vector, a named numeric, or (bad cases) list-like;
# use unlist so we never leave min()/max() seeing a list from vapply().
#' @noRd
flatten_range_bound <- function(r, j) {
  if (is.null(r)) {
    return(NA_real_)
  }
  if (length(r) < j) {
    return(NA_real_)
  }
  el <- .subset2(r, j)
  flat <- suppressWarnings(as.double(unlist(
    el,
    recursive = TRUE,
    use.names = FALSE
  )))
  if (!length(flat)) {
    return(NA_real_)
  }
  flat[[1L]]
}

#' @noRd
mixture_range_extrema <- function(r_list) {
  n <- length(r_list)
  lo <- rep(NA_real_, n)
  hi <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    ri <- r_list[[i]]
    lo[i] <- flatten_range_bound(ri, 1L)
    hi[i] <- flatten_range_bound(ri, 2L)
  }
  c(
    min(lo, na.rm = TRUE),
    max(hi, na.rm = TRUE)
  )
}

#' Mixture of distributions (`mix()` with robust support bounds)
#'
#' Same interface as [distplyr::mix()]. For non-finite components, support
#' `c(r1, r2)` passed into [distionary::distribution()] uses a small helper to
#' read bounds from `range()`, so list-like bounds do not break [base::min()] /
#' [base::max()].
#'
#' @param ... Passed to `distplyr:::pair_dots_num()` (same as `mix()`).
#' @param weights,na_action_dst,na_action_w See [distplyr::mix()].
#'
#' @returns A `dst` object (see [distplyr::mix()]).
#' @export
mix2 <- function(
  ...,
  weights = 1,
  na_action_dst = c("null", "drop", "fail"),
  na_action_w = c("null", "drop", "fail")
) {
  preprocess <- distplyr:::pair_dots_num(
    ...,
    num = weights,
    na_action_dst = na_action_dst,
    na_action_num = na_action_w
  )
  if (distionary::is_distribution(preprocess)) {
    return(preprocess)
  }
  dsts <- preprocess$dsts
  weights <- preprocess$num
  if (length(weights) == 1L) {
    return(dsts[[1L]])
  }
  probs <- weights / sum(weights, na.rm = TRUE)
  all_finite <- all(vapply(
    dsts,
    function(d) distionary::pretty_name(d) == "Finite",
    FUN.VALUE = logical(1L)
  ))
  if (all_finite) {
    x <- lapply(dsts, function(d) distionary::parameters(d)[["outcomes"]])
    p <- lapply(dsts, function(d) distionary::parameters(d)[["probs"]])
    pnew <- Map(`*`, probs, p)
    l <- distionary:::aggregate_weights(unlist(x), unlist(pnew))
    return(distionary::dst_empirical(l[["y"]], weights = l[["weight"]]))
  }
  rm("weights", "preprocess")
  if (is.list(dsts[[1]]) && !distionary::is_distribution(dsts[[1]])) {
    dsts <- lapply(dsts, \(x) x[[1]])
  }
  r <- lapply(dsts, range)
  ext <- mixture_range_extrema(r)
  r1 <- ext[[1L]]
  r2 <- ext[[2L]]
  var_type <- vapply(dsts, distionary::vtype, FUN.VALUE = character(1L))
  var_unique <- unique(var_type)
  if (length(var_unique) == 1L) {
    v <- var_unique
  } else if ("unknown" %in% var_type) {
    v <- "unknown"
  } else {
    v <- "mixed"
  }
  d <- distionary::distribution(
    cdf = function(x) {
      cdf_vals <- lapply(dsts, distionary::eval_cdf, at = x)
      p_times_cdfs <- mapply(`*`, probs, cdf_vals, SIMPLIFY = FALSE)
      Reduce(`+`, p_times_cdfs)
    },
    survival = function(x) {
      surv_vals <- lapply(dsts, distionary::eval_survival, at = x)
      p_times_survs <- mapply(`*`, probs, surv_vals, SIMPLIFY = FALSE)
      Reduce(`+`, p_times_survs)
    },
    pmf = function(x) {
      pmf_vals <- lapply(dsts, distionary::eval_pmf, at = x)
      p_times_f <- mapply(`*`, probs, pmf_vals, SIMPLIFY = FALSE)
      Reduce(`+`, p_times_f)
    },
    density = function(x) {
      density_vals <- lapply(dsts, distionary::eval_density, at = x)
      p_times_f <- mapply(`*`, probs, density_vals, SIMPLIFY = FALSE)
      Reduce(`+`, p_times_f)
    },
    realise = function(n) {
      if (n == 0L) {
        return(numeric())
      }
      k <- length(dsts)
      id <- sample(seq_len(k), size = n, replace = TRUE, prob = probs)
      vapply(
        id,
        function(i) distionary::realise(dsts[[i]], n = 1L),
        FUN.VALUE = numeric(1L)
      )
    },
    range = c(r1, r2),
    .vtype = v,
    .name = "Mixture",
    .parameters = list(distributions = dsts, probs = probs)
  )
  distionary:::new_distribution(d, class = "mixture")
}
