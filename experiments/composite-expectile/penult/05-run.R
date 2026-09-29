# ---------------------------------------------------------------------------
# Effective threshold and kernel for every estimator in this project, across
# seven populations. Population level throughout: no Monte Carlo.
#
# Validation of the pipeline against the smooth-graft project's numbers, which
# were produced independently: MLE's effective threshold there is R = 13-16 and
# its kernel spans R = 6.5 to 142. Here MLE gives 11.5-15.5 across populations,
# a kernel 10% point of 6.3-6.7 everywhere, and 141.4 for the 90% point on
# LP3(0.6, 0.2). The scale check (the target's own sigma against the penultimate
# GPD at R_e carried back to u) holds to 2.3%, against their 2.5%.
# ---------------------------------------------------------------------------
source("penult/04-setup.R")
library(parallel)
E0 <- seq(0.05, 13, length.out = 60); HW <- 0.25

## ---- the estimators, as closures of the exceedance law ---------------------
# cc and s are recomputed from the perturbed law inside the kernel, because the
# estimator itself sets them from the data before optimising; holding them fixed
# would measure a different procedure.
# Every estimator starts from the closed-form L-moment target of the law it is
# handed. A fixed numeric start is not safe across populations: on LP3(0.6,0.2)
# the exceedance scale is 0.19, and starting the Newton at sigma = 1 sends every
# composite estimator to the shape boundary at -0.84 instead of the root at 0.21.
# The formal is .kind, not kind: R partially matches named arguments to formals
# that sit before `...`, so an estimator argument named k silently bound to
# `kind` and switch() dispatched on a number -- the inverted Huber came back as
# maximum likelihood, byte for byte, with no error anywhere.
mk_fit <- function(.kind, ...) { args <- list(...); force(.kind)
  function(ex, start = target_lmom(ex)) switch(.kind,
    mle  = target_mle(ex, start),
    lmom = target_lmom(ex),
    mwle = target_mwle(ex, args$w, start),
    extr = target_extremile(ex, args$w, r_max = args$rmax, start = start),
    comp = target_composite(ex, args$loss, args$w, GRID, start,
             cc    = if (identical(args$loss, "onesided")) knot_c(ex, args$k) else NULL,
             alpha = args$alpha,
             s_el  = if (identical(args$loss, "elastile"))
                       elastile_scale(ex, args$w, args$ref) else NULL))
}

build_ests <- function(ref) {
  E <- list("MLE" = mk_fit("mle"), "L-moments" = mk_fit("lmom"))
  for (m in c(0, 2, 6, 20))
    E[[sprintf("L1 (pinball) p^%d", m)]] <- mk_fit("comp", loss = "quantile", w = pm(m))
  for (m in c(0, 1, 2, 4, 6, 10, 20))
    E[[sprintf("L2 (expectile) p^%d", m)]] <- mk_fit("comp", loss = "expectile", w = pm(m))
  for (p0 in c(0.50, 0.80, 0.90, 0.95))
    E[[sprintf("L2 smoothstep p0=%.2f", p0)]] <- mk_fit("comp", loss = "expectile", w = ssw(p0))
  E[["elastile a=0.5 p^6"]] <- mk_fit("comp", loss = "elastile", w = pm(6), alpha = 0.5, ref = ref)
  for (k in c(1, 4, 16))
    E[[sprintf("inv Huber k=%d p^6", k)]] <- mk_fit("comp", loss = "onesided", w = pm(6), k = k)
  for (m in c(2, 6)) E[[sprintf("MWLE p^%d", m)]] <- mk_fit("mwle", w = pm(m))
  for (rm in c(20, 50, 200)) E[[sprintf("extremile rmax=%d", rm)]] <- mk_fit("extr", w = pm(6), rmax = rm)
  E
}

B <- make_bases()
one_pop <- function(nm) {
  b <- B[[nm]]; ex0 <- excess_on(b)
  ref <- target_mle(ex0)                       # the elastile's scale reference
  E <- build_ests(ref)
  out <- lapply(names(E), function(en) {
    kk <- try(kernel(b, E[[en]], e0 = E0, hwid = HW), silent = TRUE)
    bad <- inherits(kk, "try-error") || any(!is.finite(kk$target))
    if (bad) return(list(row = data.frame(pop = nm, est = en, sigma = NA_real_,
        xi = NA_real_, R_e = NA_real_, sig_err = NA_real_, intK = NA_real_,
        q10 = NA_real_, q25 = NA_real_, q50 = NA_real_, q75 = NA_real_,
        q90 = NA_real_, note = "fit failed", stringsAsFactors = FALSE), kk = NULL))
    et <- try(effective_threshold(b, kk$target[2], kk$target[1]), silent = TRUE)
    if (inherits(et, "try-error")) et <- list(R = NA_real_, sigma_err = NA_real_, note = "no crossing")
    q <- try(kernel_quantiles(kk), silent = TRUE)
    if (inherits(q, "try-error")) q <- rep(NA_real_, 5)
    list(row = data.frame(pop = nm, est = en, sigma = kk$target[1], xi = kk$target[2],
           R_e = et$R, sig_err = et$sigma_err, intK = kk$total,
           q10 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q90 = q[5],
           note = et$note, stringsAsFactors = FALSE), kk = kk)
  })
  names(out) <- names(E)
  list(tab = do.call(rbind, lapply(out, `[[`, "row")),
       kern = lapply(out, `[[`, "kk"))
}

t0 <- Sys.time()
RES <- mclapply(names(B), function(nm) { r <- one_pop(nm); cat("done:", nm, "\n"); r },
                mc.cores = min(4, detectCores()))
names(RES) <- names(B)
TAB <- do.call(rbind, lapply(RES, `[[`, "tab"))
saveRDS(list(TAB = TAB, RES = RES, E0 = E0, U_RP = U_RP), "out/penult.rds")
cat(sprintf("\ntotal %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
write.csv(TAB, "out/penult.csv", row.names = FALSE)
