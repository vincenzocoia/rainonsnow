# The composite expectile estimator

A self-contained simulation study. Nothing here touches the rain-on-snow
analysis; it only borrows the repository.

## The idea

Fit a parametric family indexed by `theta` by minimising the expectile loss
integrated over levels, against a weight that is zero on the body of the
distribution and rises to one in the tail:

```
  S(theta) = sum_i  int_0^1  w(p) |p - I(y_i < E(p | theta))| (y_i - E(p | theta))^2  dp
```

where `E(. | theta)` is the family's **expectile** function. A named parametric
family is essentially always misspecified, and a global fitting criterion
(likelihood, L-moments) has to spread that misspecification across the whole
distribution. If only the tail matters, weighting the criterion towards the tail
should buy a better tail fit. The same construction with the pinball loss and
the quantile function is the composite quantile estimator from Vincenzo's
thesis; the question here is whether replacing L1 by L2 tames the variance
explosion that made the L1 version not worth using.

### A note on the loss

The expectile check function needs the **absolute value**:
`|p - I(y < e)| (y - e)^2`. Without it the criterion is unbounded below — the
factor `(p - I)` is negative for every `y` below `e`, so pushing `e` upward
drives the objective to `-Inf`. `scripts/98-validate.R` demonstrates this
numerically. The ordering `(y - E)` is right as written; it is squared, so it
does not matter, but the indicator `I(y < E)` does and is correct. The pinball
loss `(p - I(y < q))(y - q)` takes no absolute value and is already non-negative.

## The experiment

**Truth.** `F(x) = F_GEV(x; 0, 1, 0.2) * Phi((x - 1.5)/0.8)` — the cdf of
`max(GEV, Normal)` with the two independent, so a product of cdfs, as specified.
The normal sits in the middle: it dominates the body and vanishes in the tail, so
the GEV is the *exactly correct tail model* and a *wrong body model*.

**Weight.** A smoothstep, zero below `p0 = 0.95` and one above `p1 = 0.98`,
placed where the Gaussian contamination dies out (it leaves the quantiles ~0.6%
too high at `p0` and ~0% at `p1`).

**Estimators.** GEV by MLE, GEV by L-moments, composite quantile (L1),
composite expectile (L2). 2000 replicates at n = 50, 100, 200. Return periods 2
to 1000 on a log grid.

## What happened

Full write-up with figures: **[report/report.html](report/report.html)**
(published at <https://claude.ai/code/artifact/e6ef50f7-30d2-4dc8-90fb-f96889c36bbe>).
Summary:

**1. The loss needs an absolute value.** `|p - I(y<E)| (y-E)^2`. Without it the criterion is
unbounded below. `y` and `E` are the right way round; the indicator is right too.

**2. Both composite estimators remove the bias.** MLE and L-moments underestimate the 1000-year
return level by 39% and 32% asymptotically under this contamination. L1 cuts that to 0.08%, L2 to
0.7%.

**3. L2 fixes the variance explosion.** At matched weight, MSE(L2)/MSE(L1) for the 1000-year level
is 0.17-0.58 across every weight tested. L2's far-tail sd barely moves as the weight changes
(5.8-6.7); L1's swings from 4.3 to 15.0.

**4. Expectiles are anchored to the mean.** The p-expectile solves `k phi(x) + m - x = 0`, and `m`
sees the contaminated body, so the truth's expectiles differ from its GEV component's at *every*
level, decaying only like `(1-p)^xi` (still 9.7% at p = 0.98 where the quantile discrepancy is
3e-5). No weight makes L2 asymptotically unbiased; L1 is exactly unbiased.

**5. The weight was starting too late.** p0 = 0.95 was placed where the *quantile* contamination
dies out - correct for L1, badly wrong for L2. With p0 near the median the expectile version ties
the GEV MLE from T = 200 to 1000 and recovers the shape almost exactly (median xi 0.209 against a
true 0.20, where MLE gives 0.054).

**6. There is a setting where L2 wins outright.** Making the family *more* wrong does not do it:
under a graft with an exactly-GEV tail and a hostile body, the MLE is 55% low at T = 1000 and every
composite estimator does *worse*, because the heavier tail inflates their variance as fast as it
inflates the MLE's bias. What works is a near-degenerate body -
`max(GEV(0,1,0.2), N(2.5,0.3))`, i.e. tightly clustered snowmelt maxima with a heavy rain-on-snow
tail. There L2's median absolute error at T = 1000 is 0.58x the MLE's and it is closer on 84% of
datasets; L1 is *worse* than the MLE (38%). The contest is governed by
`(xi_true - xi_MLE) / sd(xi_hat_composite)`, not by the size of the bias alone.

**7. Mixing the losses helps.** A convex combination of the check functions elicits the
**alpha-elastile**, the root of
`(2a/s)[(1-2p) phi(t) + (1-p)(t-m)] + (1-a)[F(t)-p] = 0`, which always lies between the quantile
and the expectile. At a well-placed weight an interior alpha beats *both* pure losses by 5-23% on
MSE at every return period up to 500, with the best alpha rising with T (about 0.2 near the record
length, pure L2 an order of magnitude beyond). At p0 = 0.95 the effect is largely absent.

**8. Admissible weights.** For `w(p) ~ (1-p)^(-a)`: finite iff `a < 2 - xi` (L1) or `a < 2 - 2 xi`
(L2). The remembered `a < 1` is the L1 condition at xi = 1. The binding condition is on the truth:
L1 needs a finite mean (xi < 1), L2 a finite variance (xi < 1/2), above which the raw criterion is
+Inf everywhere and only the loss-difference form is usable. Weights need not be bounded by 1.

**9. Do not weight past the data.** A weight left at 1 up to p = 1 puts 27% of its mass (n = 100)
above the largest observation, where the empirical functional has saturated; the fitted shape is
dragged hard negative in proportion to that share.

Notes on `probaverse/distionary` branch `expectiles`, including a bug that returns silently wrong
expectiles for xi >~ 0.5, are in [`DISTIONARY-NOTES.md`](DISTIONARY-NOTES.md).

## Layout

```
R/gev.R          GEV cdf, quantile, mean, closed-form partial moment, expectiles
R/dgp.R          the contaminated truth, plus its quantiles and expectiles
R/weight.R       the weight function
R/quadrature.R   panelled Gauss-Legendre over levels
R/estimators.R   the four estimators
R/config.R       the settings the scripts share
src/             C++ port of the expectile solver
scripts/00-design.R             DGP and weight design, asymptotic targets
scripts/01-simulate.R           the Monte Carlo
scripts/02-figures.R            summaries and figures
scripts/03-weight-sensitivity.R how the answer moves with p0
scripts/06-loss-admissibility.R which weights are admissible, and for which tails
scripts/07-dgp-picture.R        densities and frequency-magnitude of the truth
scripts/08-tail-focus-benchmark.R  sweep the weight from w == 1 to strongly tail-focused
scripts/09-elastile.R           the alpha-elastile: population behaviour and an alpha sweep
scripts/10-alt-dgp.R            the heavier max() setting
scripts/11-graft.R              a grafted truth: exactly-GEV tail, hostile body
scripts/12-sharp-contrast.R     the near-degenerate body, where L2 wins
scripts/98-validate.R           every correctness check, re-runnable
report/build_report.py          builds the standalone HTML report
out/                            results and figures
```

Run in order from this directory; each script starts with `source("R/setup.R")`.

## Requirements

R with `Rcpp`. `distionary` (branch `expectiles`) is optional — it is used only
as an independent cross-check in `scripts/98-validate.R`.
