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

Short answer: **L2 does fix the variance problem you diagnosed — decisively — but
it is not enough to make the estimator worth using here, and it brings a new
problem of its own that L1 does not have.**

### 1. Your formula needs an absolute value

```
  |p - I(y_i < E(p|theta))| (y_i - E(p|theta))^2        <- correct
   (p - I(y_i < E(p|theta))) (y_i - E(p|theta))^2       <- as written; unbounded below
```

Without the absolute value the factor is negative for every observation below
`E`, so pushing `E` upward sends the criterion to `-Inf`; minimising it numerically
just runs to the boundary. `y` and `E` are the right way round, and the
indicator `I(y_i < E)` is right too. (The pinball loss takes no absolute value —
`(p - I)(y - q)` is already non-negative — which is probably where the slip came
from.) Demonstrated in `scripts/98-validate.R` section 1.

### 2. Expectiles are anchored to the mean, so they never fully escape the body

This is the finding I did not expect, and it is specific to L2.

The p-quantile of the contaminated truth equals the p-quantile of its GEV
component as soon as the normal's cdf reaches one — quantiles are tail-local.
Expectiles are not. The p-expectile solves `k phi(x) + m - x = 0`, and while
`phi(x) = E[(Y-x)^+]` only sees the tail, the mean `m` sees everything, including
the contaminated body. So the truth's expectiles are contaminated at *every*
level:

| level p | return period | quantile still off the GEV by | expectile still off by |
|---|---|---|---|
| 0.95 | 20 | 0.6% | 16.0% |
| 0.98 | 50 | 0.00003% | 9.7% |
| 0.99 | 100 | ~0 | 7.0% |
| 0.999 | 1000 | ~0 | 2.9% |
| 0.9999 | 10000 | ~0 | 1.4% |

The contamination leaves the quantiles at machine precision by p = 0.98 and
leaves the expectiles at roughly `(1-p)^xi` — that is, barely. So there is no
weight function that makes the composite **expectile** estimator asymptotically
unbiased under body contamination, whereas the composite **quantile** estimator
is exactly unbiased with the weight used here. Fitting to n = 2e6:

| estimator | asymptotic (mu, sigma, xi) | % bias in the 100-yr level | in the 1000-yr level |
|---|---|---|---|
| GEV MLE | (1.397, 0.890, 0.062) | -18.7% | -39.1% |
| GEV L-moments | (1.384, 0.819, 0.120) | -15.1% | -31.6% |
| Composite quantile (L1) | (0.011, 0.999, 0.200) | **+0.02%** | **-0.08%** |
| Composite expectile (L2) | (0.832, 0.879, 0.219) | +3.4% | +0.7% |

The tail-focused idea works exactly as intended — both composite estimators
essentially eliminate the 20-40% bias that MLE and L-moments carry. L1 does it
perfectly; L2 gets within a few percent, and the gap is the mean-anchoring above.

### 3. L2 roughly halves the variance explosion

This is your hypothesis, and it holds up. Standard deviation of the estimated
1000-year return level at n = 100:

| weight design | L1 | L2 |
|---|---|---|
| p0 = 0.90, uncapped | 11.05 | **6.15** |
| p0 = 0.95, uncapped | 8.73 | **5.81** |
| p0 = 0.90, capped | 26.31 | **7.33** |
| p0 = 0.95, capped | 23.81 | **7.06** |

L2's far-tail standard deviation is 40-70% below L1's in every design tested.
The MSE at T = 1000 relative to GEV MLE goes from 1.9 (L1) to 1.10 (L2) at the
primary design, and from 15.4 to 1.31 with the weight capped. That is the
variance explosion you saw, substantially tamed, exactly as the L1-to-L2
argument predicts.

### 4. It still is not worth using — but it is a much closer call

At the primary sample size (n = 100), no composite design beats the GEV MLE on
MSE anywhere from T = 20 to T = 1000. The best any of them manages in the tail is
a tie: composite expectile with `p0 = 0.98` reaches an MSE ratio of 1.05 at
T = 1000. (Two designs do dip below the MLE in a narrow band around T = 10-13 —
composite quantile with `p0 = 0.90`, ratio 0.79 — but that is inside the
contaminated region, not the tail the estimator is for. Plain L-moments also
beats the MLE for T = 3-10.) The pattern:

| n = 100, MSE relative to GEV MLE | T=20 | T=50 | T=100 | T=200 | T=500 | T=1000 |
|---|---|---|---|---|---|---|
| GEV L-moments | 1.35 | 1.30 | 1.14 | 1.10 | 1.09 | 1.11 |
| Composite quantile (L1) | 4.91 | 2.30 | 1.55 | 1.36 | 1.51 | 1.90 |
| Composite expectile (L2) | 24.02 | 9.84 | 3.58 | 2.07 | **1.28** | **1.10** |

L2 wins over L1 at the far end and loses badly in the middle — the middle is
where L2's mean-anchoring bias bites (it is biased +1.7 at T = 20 where L1 is
biased 0.0). So the two failure modes trade places: L1 is unbiased and wild, L2
is steadier and slightly biased, and they cross between T = 200 and T = 500.

The cleanest way to see the verdict is the trade each estimator makes against the
MLE at T = 1000 (n = 100), in MSE units:

| | bias^2 | variance | MSE | bias^2 bought | variance paid |
|---|---|---|---|---|---|
| GEV MLE | 33.5 | 8.2 | 41.7 | — | — |
| Composite quantile (L1) | 3.2 | 76.1 | 79.3 | -30.3 | +67.9 |
| Composite expectile (L2) | 12.1 | 33.7 | 45.8 | -21.4 | +25.5 |

L1 pays 2.2 units of variance for every unit of squared bias it removes. L2 pays
1.2. So switching to L2 turns a badly losing trade into a nearly even one — real,
and the direction you predicted — but it does not flip it.

The ordering is not permanent; it is a sample-size race. Where each estimator
beats the MLE:

| | n = 50 | n = 100 | n = 200 |
|---|---|---|---|
| GEV L-moments | T in [2.3, 9.8] | T in [2.8, 9.8] | **T in [2.8, 1000]** |
| Composite quantile (L1) | nowhere | nowhere | T in [107, 384], best 0.85 at T = 203 |
| Composite expectile (L2) | nowhere | nowhere | nowhere |

By n = 200 the composite quantile estimator does start beating the MLE in the
tail — and, worth noting, plain L-moments beats the MLE at *every* return period
by then. The composite expectile estimator never does, because by the time its
variance has come down enough its mean-anchoring bias is what is left.

### 5. Admissible weights — and one condition that hits L2 much harder than L1

You warned that the weight cannot be arbitrary. That is right, and there turn out
to be two separate conditions. `scripts/06-loss-admissibility.R` derives both and
checks them numerically.

**Condition on the weight.** Writing `u = 1 - p`, the population risk at level p
decays as

```
  quantile  (L1):   r_p ~ u^(1 - xi)
  expectile (L2):   r_p ~ u^(1 - 2 xi)
```

(measured slopes converge to these to three decimals as u -> 0). So for a weight
that blows up like `w(p) ~ (1-p)^(-a)`, the criterion is finite iff

```
  quantile  (L1):   a < 2 - xi
  expectile (L2):   a < 2 - 2 xi
```

At xi = 0.2 that is `a < 1.8` and `a < 1.6` — so a weight going like `1/(1-p)` is
still admissible for both, and `1/(1-p)^2` for neither. This is close to the
result you half-remembered, and note the L2 threshold is the tighter of the two.

The weight used here is **bounded** — identically 0 below p0 = 0.95, identically
1 above p1 = 0.98 — so `1 - w` does not decay slowly, it is exactly zero. It
clears both conditions with room to spare for every xi the estimators visit. The
weighted population risk at the true tail is 0.0048 (L1) and 0.0285 (L2), both
finite.

**Condition on the truth, and this one is the real cost of switching to L2.**
The quantile risk needs only `E|Y| < Inf`, i.e. `xi < 1`. The expectile risk
needs `E[Y^2] < Inf`, i.e. `xi < 1/2`, because `E[(e - Y)^2 1{Y < e}]` is
infinite for every `e` once the variance is:

| xi | risk at p = 0.99, quantile | expectile |
|---|---|---|
| 0.30 | 0.146 | 1.914 |
| 0.45 | 0.285 | 24.799 |
| 0.50 | 0.364 | **Inf** |
| 0.60 | 0.623 | **Inf** |

No weight function repairs this — the criterion is `+Inf` at every theta. So the
L2 version is undefined exactly in the heavy-tailed regime the method exists to
serve, which for flood frequency is not a remote corner of the parameter space.

The standard repair is to minimise the loss *difference* against a fixed
reference function,

```
  sum_i int w(p) [ rho_p(y_i, T(p | theta)) - rho_p(y_i, T_ref(p)) ] dp,
```

which has the same minimiser and is finite whenever `E|Y| < Inf`. In this
experiment xi = 0.2, so the raw criterion is finite and nothing was needed — but
it should be built in before the estimator is used on anything heavier.

### 6. A design rule that matters more than L1 vs L2

A weight that stays at one all the way to p = 1 puts a large share of its mass at
levels *above the largest observation* — 27% of it at n = 100, 45% at n = 50.
Out there the empirical quantile and expectile functions have both saturated at
the sample maximum, so the loss can only pull the fitted tail down. The result is
a fitted shape parameter dragged hard negative:

| n | share of weight above the largest observation | median fitted xi, L1 | L2 |
|---|---|---|---|
| 50 | 45% | -0.82 | -0.52 |
| 100 | 27% | -0.24 | -0.19 |
| 200 | 13% | -0.01 | -0.06 |
| 500 | 5% | 0.17 | 0.09 |
| 10 000 | 0.3% | 0.21 | 0.21 |

(True xi = 0.2.) The two columns track the first almost exactly, and both
converge to 0.2 — so this is finite-sample behaviour, not a bug. In the main run
at n = 50 this pins 89% of L1 fits and 56% of L2 fits onto the shape bound.
Capping the weight at `p = 1 - 1/n` fixes the shape (median xi goes from -0.26 to
0.19 for L1, -0.20 to 0.06 for L2) — but for L1 it then *worsens* the far-tail
MSE by an order of magnitude, because the downward drag had been acting as
accidental shrinkage. L2 survives capping intact. Whatever else you do with the
estimator, the weight should not extend past where the data can speak.

### Suggestion

If you want the tail-focus idea to actually win on MSE, the obstacle is variance,
not bias — both composite estimators already solve the bias problem. L2 buys back
about half of the variance, which is real progress but not sufficient. Two things
look more promising than further loss-function tinkering: shrinking the
effective number of free parameters over the weighted region (the loss surface is
a near-flat ridge in `(mu, sigma, xi)` — six multi-starts agree to five decimals
on a loss that varies in the fourth decimal across xi from -0.15 to +0.45), and
capping the weight support as above. It would also be worth checking whether the
mean-anchoring bias can be removed by using the *model's* expectile function
against a mean-corrected target. And if the estimator is going anywhere near
heavy tails, the loss-difference form in section 5 needs to go in first — the raw
L2 criterion is infinite for xi >= 1/2.

Notes on the state of `probaverse/distionary` (branch `expectiles`), including a
bug found while cross-checking against it, are in
[`DISTIONARY-NOTES.md`](DISTIONARY-NOTES.md).

## Figures

| file | what |
|---|---|
| `out/fig-design.png` | the truth, its GEV tail, and where the weight switches on |
| `out/fig-headline-return-level.png` | return period vs return level, median and central 95%, all six fits |
| `out/fig-headline-mse.png` | MSE and MSE ratio over return period |
| `out/fig-mse.png` | the same for the four estimators as specified |
| `out/fig-decomposition.png` | squared bias, variance, and median absolute error |
| `out/fig-samplesize.png` | MSE ratio at n = 50, 100, 200 |
| `out/fig-shape.png` | fitted shape parameter against sample size |
| `out/fig-weight-design.png` | MSE across the five weight designs |
| `out/fig-admissibility.png` | how each population risk decays as p -> 1 |

Numbers behind all of it: `out/design.txt`, `out/results.txt`,
`out/weight-design.txt`, `out/shape-diagnostic.txt`, `out/headline.txt`,
`out/admissibility.txt`, `out/validation.txt`.

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
scripts/98-validate.R           every correctness check, re-runnable
out/                            results and figures
```

Run in order from this directory; each script starts with `source("R/setup.R")`.

## Requirements

R with `Rcpp`. `distionary` (branch `expectiles`) is optional — it is used only
as an independent cross-check in `scripts/98-validate.R`.
