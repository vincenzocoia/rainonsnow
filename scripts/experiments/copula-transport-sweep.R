# Where does copula transport beat the mixture, and where does it not?
#
# copula-transport-marginal.R reports one configuration in detail. This sweeps
# the things that should govern the answer:
#
#   rho     -- dependence strength, i.e. how much of the marginal's tail
#              heaviness the predictor explains. CEVI = 1 - rho^2 for the
#              Gaussian copula, so rho = 0 means transport can add nothing and
#              rho -> 1 means it can add almost everything.
#   n       -- sample size. The mixture's tail peels away from the truth at a
#              point that depends on n, so its deficit should shrink slowly as
#              n grows while transport's should not depend on n the same way.
#   n_bins  -- how finely the covariate is sliced. More bins means a purer
#              conditional but fewer points to fit each one.
#
#   Rscript scripts/experiments/copula-transport-sweep.R

SWEEP_MODE <- TRUE
repo <- getwd()
while (!file.exists(file.path(repo, "DESCRIPTION")) && dirname(repo) != repo) {
  repo <- dirname(repo)
}
source(file.path(repo, "scripts", "experiments", "copula-transport-marginal.R"))

set.seed(4242)

N_REP <- 150L
REPORT_P <- c(1e-3, 1e-4)

flush_cat <- function(...) {
  cat(...)
  utils::flush.console()
}

# Rebuild the DGP for a given rho rather than using the file's fixed RHO.
make_dgp_copula <- function(rho, xi_y) {
  SYr <- function(y) (1 + xi_y * y)^(-1 / xi_y)
  QYr <- function(p) (p^(-xi_y) - 1) / xi_y
  function(n) {
    u <- stats::runif(n)
    v <- stats::pnorm(rho * stats::qnorm(u) + sqrt(1 - rho^2) * stats::rnorm(n))
    list(u = u, y = QYr(1 - v), true_S = SYr)
  }
}

one_config <- function(rho, n_obs, n_bins, xi_y = 0.25, n_rep = N_REP) {
  dgp <- make_dgp_copula(rho, xi_y)
  truth <- true_level(dgp(10L)$true_S, REPORT_P)
  out <- array(NA_real_, c(n_rep, length(REPORT_P), 4L))
  shp <- matrix(NA_real_, n_rep, 2L)
  for (r in seq_len(n_rep)) {
    d <- dgp(n_obs)
    cf <- fit_conditionals(d$u, d$y, n_bins = n_bins, scale_model = "free")
    cf_s <- fit_conditionals(d$u, d$y, n_bins = n_bins, scale_model = "loglinear")
    shp[r, ] <- c(cf$shape, cf_s$shape)
    out[r, , 1] <- est_direct_pot(d$u, d$y, REPORT_P) / truth
    out[r, , 2] <- est_dl_mixture(cf, REPORT_P) / truth
    out[r, , 3] <- est_transport(cf, d$u, d$y, REPORT_P) / truth
    out[r, , 4] <- est_transport(cf_s, d$u, d$y, REPORT_P) / truth
  }
  summ <- function(m, j) {
    e <- out[, j, m]
    e <- e[is.finite(e)]
    if (!length(e)) {
      return(c(NA, NA))
    }
    c(stats::median(e), sqrt(mean(log(e)^2)))
  }
  list(
    truth = truth,
    cond_shape = mean(shp[, 1], na.rm = TRUE),
    cond_shape_smooth = mean(shp[, 2], na.rm = TRUE),
    stats = lapply(seq_len(4L), function(m) {
      vapply(seq_along(REPORT_P), function(j) summ(m, j), numeric(2L))
    })
  )
}

methods <- c("direct_pot", "dl_mixture", "transport", "transport+smooth")

cat("\n=====================================================================\n")
cat("SWEEP  median ratio to truth, and RMSE of the log ratio\n")
cat(sprintf("       %d replicates, true marginal EVI 0.25\n", N_REP))
cat("       'cond EVI fitted' is the free-scale fit; transport+smooth uses the\n")
cat("       log-linear scale model instead (true value = CEVI x 0.25).\n")
cat("=====================================================================\n")

header <- function() {
  cat(sprintf(
    "\n%-6s %-6s %-6s %-8s %-12s %s\n",
    "rho",
    "n",
    "bins",
    "CEVI",
    "cond EVI",
    paste(sprintf("%-26s", paste0("p=", format(REPORT_P))), collapse = "")
  ))
  cat(sprintf(
    "%-6s %-6s %-6s %-8s %-12s %s\n",
    "",
    "",
    "",
    "",
    "fitted",
    paste(rep(sprintf("%12s%14s", "median", "RMSE(log)"), length(REPORT_P)), collapse = "")
  ))
  cat(strrep("-", 100), "\n")
}

# bins = 10 at n = 300 left too few points per bin and made everything worse,
# so the sweep stays at 5 and varies the things that matter: dependence
# strength and sample size.
grid <- expand.grid(
  rho = c(0.5, 0.7, 0.9),
  n_obs = c(300L, 1000L, 3000L),
  n_bins = 5L,
  KEEP.OUT.ATTRS = FALSE
)

header()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  res <- one_config(g$rho, g$n_obs, g$n_bins)
  for (m in seq_along(methods)) {
    st <- res$stats[[m]]
    cells <- paste(
      vapply(
        seq_along(REPORT_P),
        function(j) sprintf("%11.2fx%14.2f", st[1, j], st[2, j]),
        character(1L)
      ),
      collapse = ""
    )
    if (m == 1L) {
      cat(sprintf(
        "%-6.1f %-6d %-6d %-8.2f %-12.3f %s   %s\n",
        g$rho,
        g$n_obs,
        g$n_bins,
        1 - g$rho^2,
        res$cond_shape,
        cells,
        methods[m]
      ))
    } else {
      cat(sprintf(
        "%-6s %-6s %-6s %-8s %-12s %s   %s\n",
        "",
        "",
        "",
        "",
        "",
        cells,
        methods[m]
      ))
    }
  }
  cat("\n")
  flush_cat("")
}

cat("=====================================================================\n")
cat("Read the median column for bias and the RMSE column for total error.\n")
cat("The transport's case rests on rho: at rho = 0.3 the copula explains\n")
cat("almost none of the tail heaviness (CEVI 0.91) and there is nothing to\n")
cat("gain; at rho = 0.9 it explains most of it (CEVI 0.19).\n")
cat("=====================================================================\n")
