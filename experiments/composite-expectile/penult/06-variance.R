# ---------------------------------------------------------------------------
# Q3. The fair comparison.
#
# A composite estimator fitted above R = 4 reads the tail as if it had been
# handed a threshold at R_e. Raising a hard threshold to R_e would achieve the
# same bias -- and throw away most of the sample. So: does the composite get
# there more cheaply than the threshold does?
#
# "At the same level" admits two readings and they differ by a factor of ~3.5,
# because maximum likelihood has an effective threshold of its own. Both are
# reported:
#   (a) NOMINAL match: MLE refitted with its hard threshold at R_e.
#   (b) EFFECTIVE match: MLE refitted at whatever nominal threshold makes ITS
#       effective threshold equal R_e. This is the comparison at equal bias,
#       and the one that answers whether the composite is a cheaper route out.
# ---------------------------------------------------------------------------
source("penult/04-setup.R")
library(parallel)
NC <- min(4, detectCores())
POP <- pop_gev(0.2)
W <- list("p^0" = pm(0), "p^2" = pm(2), "p^6" = pm(6), "p^20" = pm(20))
RS <- c(4, 6, 8, 11, 14, 20, 30, 50, 80, 130, 200)      # nominal MLE thresholds

## ---- population: MLE's effective threshold against its nominal one ---------
mle_curve <- do.call(rbind, lapply(RS, function(R) {
  b <- make_base(POP, POP$Qs(1/R)); ex <- excess_on(b)
  th <- target_mle(ex); et <- effective_threshold(b, th[2], th[1])
  data.frame(R_nom = R, xi = th[2], R_eff = et$R)
}))
cat("MLE: nominal threshold -> effective threshold\n"); print(round(mle_curve, 3))
inv_eff <- splinefun(log(mle_curve$R_eff), log(mle_curve$R_nom), method = "monoH.FC")
nominal_for_effective <- function(Re) exp(inv_eff(log(Re)))

## ---- population: each composite estimator's own effective threshold --------
b4 <- make_base(POP, POP$Qs(1/4)); ex4 <- excess_on(b4)
CE <- list()
for (nm in names(W)) CE[[paste("L2", nm)]] <-
  list(fit = function(z, w) fit_gpd2_composite(z, GRID, w, "expectile"), w = W[[nm]],
       tgt = target_composite(ex4, "expectile", W[[nm]], GRID, target_lmom(ex4)))
CE[["invHuber k=4 p^6"]] <-
  list(fit = function(z, w) fit_gpd2_onesided(z, GRID, w, 4), w = W[["p^6"]],
       tgt = target_composite(ex4, "onesided", W[["p^6"]], GRID, target_lmom(ex4),
                       cc = knot_c(ex4, 4)))
for (nm in names(CE)) {
  et <- effective_threshold(b4, CE[[nm]]$tgt[2], CE[[nm]]$tgt[1])
  CE[[nm]]$R_e <- et$R
}

## ---- Monte Carlo ----------------------------------------------------------
one <- function(y) {
  u <- as.numeric(quantile(y, 0.75, type = 1)); z <- y[y > u] - u
  o <- c(vapply(names(CE), function(nm) {
           p <- try(CE[[nm]]$fit(z, CE[[nm]]$w), silent = TRUE)
           if (inherits(p, "try-error") || anyNA(p)) NA_real_ else p[2] }, numeric(1)),
         vapply(RS, function(R) {
           uu <- as.numeric(quantile(y, 1 - 1/R, type = 1)); zz <- y[y > uu] - uu
           p <- try(gpd2_mle(zz), silent = TRUE)
           if (inherits(p, "try-error") || anyNA(p)) NA_real_ else p[2] }, numeric(1)))
  names(o) <- c(names(CE), paste0("MLE R=", RS)); o
}

run_n <- function(n, nrep, seed) {
  set.seed(seed)
  dat <- lapply(seq_len(nrep), function(i) POP$Qs(runif(n)))
  M <- do.call(rbind, mclapply(dat, one, mc.cores = NC))
  data.frame(est = colnames(M), n = n,
             mean = colMeans(M, na.rm = TRUE),
             sd = apply(M, 2, sd, na.rm = TRUE),
             rmse = sqrt(colMeans(sweep(M, 2, 0.2)^2, na.rm = TRUE)),
             nfail = colSums(is.na(M)), row.names = NULL)
}

sink("out/penult-variance.txt", split = TRUE)
cat("=== Q3: variance at matched effective threshold, GEV(0,1,0.2) ===\n\n")
cat("MLE's own effective threshold against its nominal one:\n")
print(round(mle_curve, 2))
cat("\nComposite estimators fitted above R = 4, effective thresholds:\n")
for (nm in names(CE)) cat(sprintf("  %-18s xi* = %.4f   R_e = %6.1f   (MLE would need a nominal threshold of R = %.1f)\n",
  nm, CE[[nm]]$tgt[2], CE[[nm]]$R_e, nominal_for_effective(CE[[nm]]$R_e)))

RES <- list()
for (cfg in list(c(1000, 1200, 811), c(5000, 600, 812))) {
  n <- cfg[1]; nrep <- cfg[2]
  cat(sprintf("\n\n--- n = %d, %d replicates ---\n", n, nrep))
  t0 <- Sys.time(); D <- run_n(n, nrep, cfg[3]); RES[[as.character(n)]] <- D
  cat(sprintf("(%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  print(format(D[, c("est","mean","sd","rmse","nfail")], digits = 4))
  ## where does MLE's sd match each composite's?
  mi <- match(paste0("MLE R=", RS), D$est)
  sdm <- D$sd[mi]; rmm <- D$rmse[mi]
  interp <- function(v, from) { o <- order(from)
    if (!isTRUE(v >= min(from, na.rm = TRUE) && v <= max(from, na.rm = TRUE))) return(NA_real_)
    exp(approx(log(from[o]), log(RS[o]), log(v))$y) }
  cat("\n  MLE's sd(xi-hat) by nominal threshold:\n    ")
  cat(paste(sprintf("R=%g:%.4f", RS, sdm), collapse = "  "), "\n")
  cat("\n  estimator          sd(xi)   R_e   equal-bias MLE   MLE nominal R matching this sd\n")
  for (nm in names(CE)) {
    sv <- D$sd[D$est == nm]
    cat(sprintf("  %-18s %.4f %6.1f %10.1f %18s\n", nm, sv, CE[[nm]]$R_e,
        nominal_for_effective(CE[[nm]]$R_e),
        ifelse(is.na(interp(sv, sdm)), "(off the grid)", sprintf("%.1f", interp(sv, sdm)))))
  }
  cat("\n  Read: if the last column exceeds the equal-bias column, the composite reaches\n")
  cat("  that effective threshold for LESS variance than raising a hard threshold would.\n")
}
sink()
saveRDS(list(RES = RES, CE = CE, mle_curve = mle_curve, RS = RS), "out/penult-variance.rds")
cat("\nwrote out/penult-variance.txt\n")
