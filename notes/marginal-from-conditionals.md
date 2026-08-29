# Recovering a marginal tail from conditional ones

*Working note, August 2026. Code and reproducible experiments referenced
throughout; every number below is printed by a script in
`scripts/experiments/` with its output kept in `scripts/experiments/results/`.*

---

## 1. The problem: a mixture of learned conditionals cannot reach the marginal

The pipeline learns one predictive distribution of runoff per peak hour
(a quantile regression forest per grid cell) and forms the cell marginal as the
equal-weight mixture of them. This is the law of total probability over the
covariate values the sample happens to contain, and it is biased low in a way
that does not go away.

Write $U = F_X(X)$, $V = F_Y(Y)$, and let $\mathrm{CEVI}(u)$ be the extreme
value index of $1/(1-V)$ given $U = u$ — a property of the copula alone, so that
$\mathrm{EVI}(Y \mid X = x) = \mathrm{CEVI}(u)\cdot\mathrm{EVI}(Y)$. Every
component of the mixture carries the *conditional* index, so a finite mixture of
them is regularly varying with that index too: lighter than the truth by exactly
$\mathrm{CEVI}$ whenever $\mathrm{CEVI} < 1$. The marginal only acquires its full
heaviness from covariate values beyond any sample, and no amount of learning
reaches them.

Median ratio to the true $10^{-4}$ return level, Gaussian copula, generalized
Pareto marginal (`copula-transport-sweep.R`):

| $\rho$ | CEVI | $n = 300$ | $n = 1000$ | $n = 3000$ |
|---|---|---|---|---|
| 0.5 | 0.75 | 0.84x | 0.88x | 0.89x |
| 0.7 | 0.51 | 0.67x | 0.75x | 0.78x |
| 0.9 | 0.19 | 0.55x | 0.59x | 0.61x |

Tripling the sample barely moves it. The deficit is structural, not estimation
error.

## 2. What fixes it, and how the pieces are combined

With $h(v \mid u) = \partial C(u,v)/\partial u$, the conditional and the marginal
are related by $F_{Y|X=x}(y) = h\!\left(F_Y(y) \mid u\right)$. Inverting $h$ in
its first argument turns any *one* conditional into an estimate of the whole
marginal:

$$F_Y(y) = h^{-1}\!\left(F_{Y|X=x}(y) \mid u\right), \qquad u = F_X(x).$$

Two points about how this is done here.

**It is the distribution that is transported, not a tail index.** Estimating
$\mathrm{EVI}(Y)$ as $\hat\xi_{\text{cond}}/\mathrm{CEVI}$ discards the
distribution and keeps one number. Applying $h^{-1}$ to the whole conditional
survival curve keeps everything, and yields one estimate of the entire marginal
*per covariate value* — the covariate is carried through rather than integrated
out.

**Combining those estimates is an averaging problem, not a mixing one.** That
distinction decides which operations are legal. Mixing conditionals is the law
of total probability and admits only the weighted mean. Once each curve is
already an estimate of the *same* marginal, the pointwise median becomes
available, and it is well defined as a distribution: a pointwise median of
monotone survival curves is itself a monotone survival curve. It matters. Ratio
to truth, $n = 3000$, with a tail shape fitted separately at each covariate value
(`median-of-marginals.R`, part 3):

| combination | $10^{-3}$ | $10^{-4}$ | $10^{-5}$ |
|---|---|---|---|
| pointwise median | 0.97x | 0.93x | 0.90x |
| mean | 1.08x | 1.39x | **2.10x** |
| weighted mean | 1.07x | 1.39x | 2.07x |

The mean is dragged upward by whichever covariate value drew the heaviest
shape — the same domination that makes an unpooled mixture inherit
$\max_i \xi_i$, reappearing one level up. At $n = 1000$ the mean reaches 3.80x.

All arithmetic in the far tail is done in survival space. Composing $h$ with
$h^{-1}$ through probability space forms $1 - s$ for a tiny $s$, which rounds to
one; at $s = 10^{-10}$ the round trip returns a value 130 times too large.

## 3. The process these claims are measured on

$F_X$ and $Y \mid X$ are chosen and the marginal is **derived** from them, so no
part of the marginal is picked for convenience:

$$X \sim \text{Pareto}(4), \qquad Y \mid X = x \;\sim\; \mathrm{GP}(\text{scale}=x,\ \xi = 0.1),$$

that is $Y = X\cdot W$ with $W$ generalized Pareto and independent of $X$. Every
conditional is *exactly* generalized Pareto. Breiman's lemma gives the marginal
the covariate's index, $1/4$, so $\mathrm{CEVI} = 0.4$. The derived marginal is
generalized Pareto at no level a sample reaches:

| $P(Y>y)$ | $10^{-1}$ | $10^{-2}$ | $10^{-3}$ | $10^{-4}$ | $10^{-5}$ |
|---|---|---|---|---|---|
| local index $-\mathrm{d}\log S/\mathrm{d}\log y$ | 0.518 | 0.317 | 0.266 | 0.253 | 0.250 |

A generalized Pareto fitted in the top few percent reads a slope that has not
settled and then extrapolates it as if it had. That error is in the *shape of
the curve*, not in how precisely it was measured, and more data does not remove
it. This is the decomposition argument in arithmetic: the part is a clean
parametric object where the whole is not.

## 4. What is **not** established

The head-to-head against the obvious alternative — one generalized Pareto fitted
to $Y$ alone, ignoring the covariate — is currently a tie where it matters.
Ratio to truth on the process above:

| | $10^{-3}$ | $10^{-4}$ | $10^{-5}$ |
|---|---|---|---|
| transport, median, $n = 1000$ | 0.91x | 0.84x | 0.74x |
| GP fitted to $Y$ alone, $n = 1000$ | 0.91x | 0.81x | 0.69x |
| transport, median, $n = 3000$ | 0.97x | 0.94x | 0.89x |
| GP fitted to $Y$ alone, $n = 3000$ | 0.99x | 0.94x | 0.85x |
| transport, median, $n = 10000$ | 1.00x | 0.99x | 0.99x |
| GP fitted to $Y$ alone, $n = 10000$ | 0.99x | 0.96x | 0.90x |

The separation only becomes clear at $n = 10{,}000$, and this is a process built
to favour the decomposition. Three reasons to be careful:

**The transport relocates an extrapolation rather than removing it.** It trades
extrapolating the marginal tail under a generalized Pareto assumption for
extrapolating a conditional tail (shallower, better determined) *plus* the
copula's corner under a family assumption. Body data discriminates between
copula families in the corner about as poorly as it discriminates between
tail families.

**The advantage and the ill-conditioning grow together.** Where the covariate
explains much of the tail heaviness (small CEVI), the conditional tail is nearly
exponential and carries little information about its own shape — at $\rho = 0.9$,
$\mathrm{CEVI} = 0.19$, every method tested failed. Where the conditional tail is
well determined, the covariate was not explaining much. There is a middle band,
and its boundaries are not yet mapped.

**Accuracy is governed by the conditional model, not by the transport.** Holding
the conditionals exact and shifting only the shape by $+0.02$ costs a factor 1.25
at $10^{-5}$ when the extrapolation spans 8 log units and 1.37 when it spans 13 —
that amplification is what CEVI measures. But estimation error swamps it: a
kernel neighbourhood is itself a small scale mixture and comes out too heavy, and
widening it from 0.015 to 0.12 on the rank scale moves the fitted conditional
shape from 19% low to 15% high and the $10^{-5}$ level from 0.84x to 1.16x.

The decomposition does not add information about the far tail. Thirty years of
data is thirty years of data. What it can do is put the parametric assumption
where it is more nearly true, and make it checkable.

## 5. What the decomposition buys that a marginal fit cannot

1. **The assumption goes where it holds.** Demonstrated above: the conditional is
   exactly generalized Pareto, the marginal is not.

2. **A built-in consistency check.** Every covariate value estimates the *same*
   marginal, so disagreement between them is a detectable internal
   contradiction. A generalized Pareto fitted to the marginal has no such thing —
   it cannot tell you it is wrong. A deliberately misspecified copula shows a
   between-value spread of 5.77 against 1.32 for the correct one.

3. **It answers the question actually being asked.** For compound extremes the
   question is not only what the 200-year runoff is, but which combination of
   rainfall and snowmelt produces it. The marginal return level is close to a
   by-product.

## 6. Open questions

- Which comparison is the one worth winning: against the mixture (settled), or
  against a direct marginal fit (not settled)?
- Is there reason to trust a copula family in its corner more than an
  extreme-value family in its tail? If not, is claim 2 above the real
  contribution rather than a better point estimate?
- If the true copula has non-constant CEVI, transporting under a fitted Gaussian
  (constant CEVI) should make the per-covariate disagreement *systematic in $u$*
  rather than random. That would turn the spread diagnostic into an empirical
  test for non-constant CEVI — worth pursuing?
- The pipeline uses two predictors, so there is no scalar $F_X(x)$; the rank of
  the forest's own predictive summary stands in for it and the fitted copula
  absorbs what that loses. Is there a better one-dimensional reduction?

---

### Where the code is

| | |
|---|---|
| transport, ensembles, combination rules | `R/copula_transport.R` |
| the derived-marginal process | `R/transport_lab.R` |
| pipeline integration (off by default) | `R/dl_transport_marginal.R`, `scripts/5-runoff_marginals.r` |
| experiments | `scripts/experiments/{median-of-marginals,copula-transport-sweep,copula-transport-by-x}.R` |
| interactive version of the experiments | `apps/8-copula-transport-lab` |
