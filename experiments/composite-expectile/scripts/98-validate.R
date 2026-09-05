# ---------------------------------------------------------------------------
# Validation of everything the experiment rests on. Run with:
#   Rscript scripts/98-validate.R
# Output: out/validation.txt
#
# distionary's expectile branch is used as an independent cross-check where it
# is installed; the checks that do not need it still run if it is absent.
# ---------------------------------------------------------------------------
source("R/setup.R")
source("R/config.R")
dir.create("out", showWarnings = FALSE)
sink("out/validation.txt", split = TRUE)
has_distionary <- requireNamespace("distionary", quietly = TRUE)

cat("=== 1. The loss functions target what they claim ===\n")
cat("Minimising the sample loss at a single level p should return that level's\n")
cat("expectile (L2) or quantile (L1). Note the absolute value in the expectile\n")
cat("loss: without it the criterion is unbounded below.\n\n")
set.seed(11); y <- qgev(runif(4e5), 0, 1, 0.2)
for (p in c(0.3, 0.7, 0.95)) {
  e_hat <- optimize(function(x) mean(abs(p - (y < x)) * (y - x)^2), c(-3, 40), tol = 1e-9)$minimum
  e_bad <- optimize(function(x) mean((p - (y < x)) * (y - x)^2), c(-3, 1e4), tol = 1e-9)$minimum
  q_hat <- optimize(function(x) mean((p - (y < x)) * (y - x)), c(-3, 40), tol = 1e-9)$minimum
  cat(sprintf("p=%.2f  argmin |p-I|(y-e)^2 = %8.4f  (true expectile %8.4f)\n", p, e_hat, gev_expectile(p,0,1,0.2)))
  cat(sprintf("        argmin (p-I)(y-e)^2 = %8.4f  <- no absolute value: runs to the boundary\n", e_bad))
  cat(sprintf("        argmin (p-I)(y-q)   = %8.4f  (true quantile  %8.4f)\n", q_hat, qgev(p,0,1,0.2)))
}

cat("\n=== 2. Analytic partial moment phi(x) = E[(X-x)^+] ===\n")
cat("against  phi(x) = int_{F(x)}^1 (Q(u) - x) du, computed numerically.\n\n")
phi_ref <- function(x, mu, sigma, xi) sapply(x, function(xx) {
  Fx <- pgev(xx, mu, sigma, xi); if (Fx >= 1) return(0)
  integrate(function(u) qgev(u,mu,sigma,xi) - xx, Fx, 1, rel.tol=1e-12, subdivisions=2000)$value })
for (par in list(c(0,1,0.2), c(0,1,0.4), c(2,3,0.1), c(0,1,-0.2))) {
  xs <- qgev(c(0.05,0.3,0.5,0.8,0.95,0.99), par[1],par[2],par[3])
  a <- gev_partial_moment(xs, par[1],par[2],par[3]); b <- phi_ref(xs, par[1],par[2],par[3])
  cat(sprintf("GEV(%.0f,%.0f,%+.2f)  max rel err %.2e\n", par[1],par[2],par[3], max(abs(a-b)/abs(b))))
}
cat("\nand at the lower endpoint, where phi must equal (mean - endpoint):\n")
for (par in list(c(0,1,0.2), c(0,1,0.4), c(0,1,0.7))) {
  lo <- gev_endpoints(par[1],par[2],par[3])[1]
  cat(sprintf("GEV(%.0f,%.0f,%+.2f)  phi(lo) = %.10f   mean - lo = %.10f\n", par[1],par[2],par[3],
              gev_partial_moment(lo,par[1],par[2],par[3]), gev_mean(par[1],par[2],par[3]) - lo))
}

cat("\n=== 3. Expectile solver ===\n")
p <- c(0.001,0.01,0.1,0.3,0.5,0.7,0.9,0.95,0.99,0.999,0.9999,0.99999)
cat("(a) identification equation  p E[(X-e)^+] = (1-p) E[(e-X)^+]  at the solution:\n")
for (par in list(c(0,1,0.2),c(0,1,0.4),c(2,3,0.1),c(0,1,-0.2),c(-1,2,0.3),c(0,1,0.7))) {
  e <- gev_expectile(p, par[1],par[2],par[3]); m <- gev_mean(par[1],par[2],par[3])
  up <- gev_partial_moment(e, par[1],par[2],par[3]); lo <- up - (m - e)
  cat(sprintf("GEV(%+.0f,%.0f,%+.2f)  max residual %.2e\n", par[1],par[2],par[3], max(abs(p*up-(1-p)*lo))))
}
cat("\n(b) C++ solver vs the R reference implementation:\n")
for (xi in c(-0.3,-0.1,0.05,0.2,0.4,0.7,0.85))
  cat(sprintf("xi=%+.2f  max rel diff %.2e\n", xi,
      max(abs(gev_expectile_std_cpp(p,xi) - gev_expectile_std(p,xi)) / abs(gev_expectile_std(p,xi)))))
cat("\n(c) location-scale equivariance, e(p; mu, sigma) == mu + sigma e(p; 0, 1):\n")
cat("max abs diff:", max(abs(gev_expectile(p,3,5,0.25) - (3 + 5*gev_expectile(p,0,1,0.25)))), "\n")
if (has_distionary) {
  cat("\n(d) vs distionary::eval_expectile (branch 'expectiles'):\n")
  for (par in list(c(0,1,0.2),c(0,1,0.4),c(2,3,0.1),c(0,1,-0.2),c(-1,2,0.3))) {
    th <- try(distionary::eval_expectile(distionary::dst_gev(par[1],par[2],par[3]), at=p), silent=TRUE)
    cat(sprintf("GEV(%+.0f,%.0f,%+.2f)  %s\n", par[1],par[2],par[3],
      if (inherits(th,"try-error")) "distionary errored"
      else sprintf("max rel diff %.2e", max(abs(gev_expectile(p,par[1],par[2],par[3])-th)/abs(th)))))
  }
  cat("\n(e) distionary on a very heavy tail (xi = 0.7), where its numerical\n")
  cat("    integration of the survival function is expected to struggle:\n")
  th <- suppressWarnings(try(distionary::eval_expectile(distionary::dst_gev(0,1,0.7), at=p), silent=TRUE))
  cat("    distionary:", if (inherits(th,"try-error")) "error" else paste(signif(th[1:6],5), collapse=" "), "\n")
  cat("    here      :", paste(signif(gev_expectile(p,0,1,0.7)[1:6],5), collapse=" "), "\n")
} else cat("\n(d) distionary not installed; skipped.\n")

cat("\n=== 4. Fast composite loss == the direct n-by-K definition ===\n")
ref_obj <- function(y, grid, w_fun, type) {
  p <- grid$p; wts <- grid$w_quad*w_fun(p); k <- wts>0; p<-p[k]; wts<-wts[k]
  P <- matrix(p, nrow=length(y), ncol=length(p), byrow=TRUE)
  function(par){ tv <- if(type=="quantile") qgev(p,par[1],exp(par[2]),par[3])
                       else gev_expectile(p,par[1],exp(par[2]),par[3])
    d <- outer(y,tv,"-"); b <- d<0
    if(type=="quantile") sum(wts*colSums((P-b)*d)) else sum(wts*colSums(abs(P-b)*d^2)) } }
set.seed(7)
for (n in c(50,100,500)) { yy <- r_true(n)
  for (ty in c("quantile","expectile")) {
    f1 <- make_composite_objective(yy,GRID,W_FUN,ty); f2 <- ref_obj(yy,GRID,W_FUN,ty)
    d <- max(sapply(1:40, function(i){ th <- c(rnorm(1,1,1), log(runif(1,.5,3)), runif(1,-.2,.6))
      abs(f1(th)-f2(th))/max(1e-12,abs(f2(th))) }))
    cat(sprintf("n=%4d %-9s max rel diff %.2e\n", n, ty, d)) } }

cat("\n=== 5. Quadrature: fitted 100-year return level vs rule refinement ===\n")
set.seed(9); samples <- lapply(1:5, function(i) r_true(100))
for (ty in c("quantile","expectile")) {
  cfg <- list(c(6,8),c(12,8),c(24,8),c(48,8),c(96,8))
  tab <- sapply(cfg, function(cc){ g <- make_level_grid(P0, cc[1], cc[2])
    sapply(samples, function(y){ f <- fit_composite(y,g,W_FUN,ty); qgev(0.99,f[1],f[2],f[3]) }) })
  colnames(tab) <- sapply(cfg, function(cc) sprintf("%dx%d", cc[1], cc[2]))
  cat(sprintf("\n%s: max relative change vs the finest rule\n", ty))
  print(signif(apply(abs(tab-tab[,ncol(tab)])/tab[,ncol(tab)], 2, max), 3))
}
cat(sprintf("\nThe experiment uses %d x %d = %d levels.\n", N_PANEL, N_GL, GRID$n_nodes))

cat("\n=== 6. Optimiser reliability: multi-start on one sample ===\n")
set.seed(4); y <- r_true(100)
for (ty in c("quantile","expectile")) {
  obj <- make_composite_objective(y, GRID, W_FUN, ty)
  st <- list(c(0,0,0.2), c(2,log(2),0.05), c(-3,log(3),0.4), c(1,log(0.7),-0.2), c(5,0,0.3))
  r <- t(sapply(st, function(s){ o <- optim(s,obj,method="Nelder-Mead",control=list(maxit=5000,reltol=1e-14))
    o <- optim(o$par,obj,method="Nelder-Mead",control=list(maxit=5000,reltol=1e-14))
    c(o$par[1], exp(o$par[2]), o$par[3], o$value) }))
  cat(sprintf("\n%s: spread of the 5 solutions (max - min)\n", ty))
  cat(sprintf("  mu %.2e  sigma %.2e  xi %.2e  loss %.2e\n",
              diff(range(r[,1])), diff(range(r[,2])), diff(range(r[,3])), diff(range(r[,4]))))
}

cat("\n=== 7. The grafted truth ===\n")
TRG <- make_graft(0, 1, 0.4, 0.90, 3, "tight")
h <- 1e-7
cat(sprintf("join u = %.4f at F = %.2f; body N(%.4f, %.4f) truncated above u\n",
            TRG$u, TRG$Fu, TRG$a, TRG$b))
cat(sprintf("density continuous at the join: f(u-) = %.8f, f(u+) = %.8f\n",
            d_graft(TRG$u - h, TRG), d_graft(TRG$u + h, TRG)))
xt <- seq(TRG$u, 60, length.out = 200)
cat(sprintf("tail is EXACTLY the GEV above the join: max |F_graft - F_GEV| = %.3e\n",
            max(abs(p_graft(xt, TRG) - pgev(xt, 0, 1, 0.4)))))
cat(sprintf("density integrates to: %.8f\n",
            integrate(function(x) d_graft(x, TRG), -Inf, TRG$u)$value +
            integrate(function(x) d_graft(x, TRG), TRG$u, Inf, rel.tol = 1e-10)$value))
ppg <- c(0.001, 0.05, 0.3, 0.6, 0.89, 0.9, 0.91, 0.99, 0.999)
cat(sprintf("q_graft inverts p_graft: max |F(Q(p)) - p| = %.3e\n",
            max(abs(p_graft(q_graft(ppg, TRG), TRG) - ppg))))
set.seed(2); yg <- r_graft(2e6, TRG)
cat(sprintf("simulated vs analytic quantiles, max rel diff: %.4f\n",
            max(abs(quantile(yg, ppg) - q_graft(ppg, TRG)) / abs(q_graft(ppg, TRG)))))
cat(sprintf("fraction above the join: simulated %.4f, analytic %.4f\n",
            mean(yg >= TRG$u), 1 - TRG$Fu))

cat("\n=== 8. The alpha-elastile ===\n")
pe <- c(0.1, 0.5, 0.9, 0.99, 0.999); sc <- 1.5
qq <- qgev(pe, 0, 1, 0.2); ee <- gev_expectile(pe, 0, 1, 0.2); mm <- gev_mean(0, 1, 0.2)
cat("limits: alpha -> 0 is the quantile, alpha -> 1 the expectile\n")
for (a in c(1e-6, 1 - 1e-6))
  cat(sprintf("  alpha = %-9g max|.-q| = %.2e  max|.-e| = %.2e\n", a,
      max(abs(gev_elastile(pe,0,1,0.2,a,sc) - qq)),
      max(abs(gev_elastile(pe,0,1,0.2,a,sc) - ee))))
cat("first-order condition G(t) = 0 at the solution:\n")
for (a in c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  tt <- gev_elastile(pe, 0, 1, 0.2, a, sc)
  G <- (2*a/sc)*((1-2*pe)*gev_partial_moment(tt,0,1,0.2) + (1-pe)*(tt-mm)) +
       (1-a)*(pgev(tt,0,1,0.2) - pe)
  cat(sprintf("  alpha = %.1f  max|G| = %.2e\n", a, max(abs(G))))
}
M <- sapply(seq(0,1,by=0.1), function(a) gev_elastile(pe,0,1,0.2,a,sc))
cat(sprintf("always between the quantile and the expectile: %s\n",
            all(M >= pmin(qq,ee) - 1e-9 & M <= pmax(qq,ee) + 1e-9)))
cat(sprintf("monotone in alpha at every level: %s\n",
            all(apply(M, 1, function(r) all(diff(r) >= -1e-9) || all(diff(r) <= 1e-9)))))

cat("\n=== 9. Truth: simulated vs analytic ===\n")
set.seed(3); ysim <- r_true(2e6)
pp <- c(0.5,0.9,0.99,0.999)
cat("quantiles  empirical:", paste(round(quantile(ysim, pp),4), collapse=" "), "\n")
cat("           analytic :", paste(round(q_true(pp),4), collapse=" "), "\n")
cat("mean       empirical:", round(mean(ysim),4), "  analytic:", round(mean_true(),4), "\n")
sink()
