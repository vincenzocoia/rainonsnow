# Shared configuration for the penultimate study.
source("R/setup.R"); source("R/config.R"); source("R/gpd.R"); source("R/gpd_estimators.R")
source("R/gpd2.R"); source("R/invhuber.R"); source("R/extremile.R")
source("penult/00-pop.R"); source("penult/01-excess.R")
source("penult/02-targets.R"); source("penult/03-kernel.R")

GRID <- make_level_grid(0, 24, 8)
U_RP <- 4                                     # fit above the population 1/4 quantile
POPS <- c(lapply(c(0, 0.2, 0.45), pop_gev),
          list(pop_lp3(0.6, 0.2), pop_lp3(1.6, 0.2),
               pop_lp3(0.6, 0.45), pop_lp3(1.6, 0.45)))
names(POPS) <- vapply(POPS, `[[`, character(1), "name")

pm <- function(m) { force(m); function(p) p^m }             # w(p) = p^m
ssw <- function(p0) make_weight(p0, 1)                      # smoothstep from p0

# The elastile's scale s balances the two loss components; take it as the ratio
# of the weighted L2 to the weighted L1 loss at a reference GPD, as the
# finite-sample study did, so alpha means the same thing across populations.
elastile_scale <- function(ex, w_fun, ref) {
  phiL <- mk_phiL(ex); p <- GRID$p; wq <- GRID$w_quad * w_fun(p)
  k <- wq > 1e-14; p <- p[k]; wq <- wq[k]
  te <- gpd_expectile(p, 0, ref[1], ref[2]); tq <- qgpd(p, 0, ref[1], ref[2])
  sum(wq * (p * ex$phic(te) + (1-p) * phiL(te))) /
  sum(wq * (p * ex$phic(tq) + (1-p) * phiL(tq)))
}
# The knot the fitter would pick: k times the IQR of the exceedances themselves.
knot_c <- function(ex, k = 4) k * (ex$Qc(0.75) - ex$Qc(0.25))

# One base grid per population, built once and reused by every estimator and
# every kernel perturbation.
make_bases <- function(Emax = 20, N = 6001)
  lapply(POPS, function(P) make_base(P, P$Qs(1/U_RP), Emax = Emax, N = N))
