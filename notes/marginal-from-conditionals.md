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

**But read the return period before reading the table.** Those figures are at
$p = 10^{-4}$ — a 10,000-year event. At the levels hydrology actually asks
about, the 20- to 200-year event, the same comparison on the process of
section 8 puts the two within about ten percent of each other; the truncation
only takes over past a factor of roughly $10^5$ in survival probability. The
structural bias is real and the mechanism is exactly as described, but it is
not what is wrong with a design flood estimate.

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

## 4. The head-to-head, and the failure it exposes

The comparison the whole argument stands on is not against the mixture — that
one is settled — but against the obvious alternative: ignore the covariate and
fit one generalized Pareto to $Y$. `transport-vs-marginal-fit.R` runs it on the
same process with $X \sim \text{Pareto}(4)$, sweeping CEVI from 0.1 to 0.8 by
moving the conditional shape while the marginal's asymptotic index stays at
0.25. Six estimators on the same sample, paired; the kernel bandwidth follows
one fixed rule and is never tuned per setting.

Median ratio to the true $10^{-5}$ level, $n = 10{,}000$:

| CEVI | 0.1 | 0.2 | 0.4 | 0.6 | 0.8 |
|---|---|---|---|---|---|
| transport, **true** copula | 1.00 | 1.00 | 1.00 | 0.98 | 0.99 |
| transport, fitted Gaussian | 0.40 | 0.43 | 0.53 | 0.61 | 0.71 |
| transport, fitted survival Clayton | 0.47 | 0.53 | 0.66 | 0.81 | 0.97 |
| mixture over covariate values | 0.59 | 0.61 | 0.67 | 0.74 | 0.80 |
| GP on $Y$ alone, 90% threshold | 0.83 | 0.79 | 0.91 | 0.94 | 1.02 |

Three readings, none of them comfortable.

**The machinery is right and the model risk is the whole story.** Given the true
copula the transport is exact — 1.00 bias, rms log error 0.14 to 0.21. Given a
copula fitted from the sample by Kendall's tau it is the *worst* estimator in
the table: worse than the mixture it was designed to replace, and worse than
ignoring the covariate entirely. Everything between those two rows is copula
misspecification.

**It is model error, not estimation error.** The fitted-Gaussian bias at
$n = 10{,}000$ is 0.40; at $n = 1000$ it is 0.39. More data does not touch it.

**The direction inverts the hoped-for story.** The transport gets *worse* as
CEVI falls. The more of the tail heaviness the covariate could explain, the more
of the answer is supplied by the corner of the copula — and Kendall's tau is a
statistic of the body.

*(A correction to how this was presented earlier: the transport rows in
`median-of-marginals.R` all use the process's own copula. They are the
upper-bound row above, not an achievable method.)*

## 5. The repair: choose the copula by making the covariate values agree

The transport supplies a criterion a marginal fit does not have. Every covariate
value transports to an estimate of the **same** marginal, so under the right
copula they must agree, and disagreement is computable *without the truth*. Take
the standard deviation of $\log S$ across covariate values, summed over three
reference levels read off the sample itself, and minimise it over a grid of
families and parameters. This targets the corner directly instead of inheriting
it from a body statistic. Results (`copula-by-agreement.R`, 60 replicates):

**It detects the misspecification.** Median disagreement at $n = 10{,}000$:

| CEVI | 0.1 | 0.2 | 0.4 | 0.6 | 0.8 |
|---|---|---|---|---|---|
| under the true copula | 0.65 | 0.65 | 0.68 | 0.66 | 0.64 |
| under the Kendall tau copula | 4.07 | 3.10 | 2.25 | 1.69 | 1.34 |

and, replicate by replicate, Spearman correlation between log disagreement and
$\lvert\log(\text{error})\rvert$ runs 0.53 to 0.74. The ensemble knows when it
is wrong.

**And it repairs the estimate — up to a point.** Root mean square log ratio at
$10^{-5}$, $n = 10{,}000$:

| CEVI | 0.1 | 0.2 | 0.4 | 0.6 | 0.8 |
|---|---|---|---|---|---|
| transport, true copula (bound) | 0.19 | 0.21 | 0.14 | 0.15 | 0.14 |
| transport, copula from Kendall tau | 0.92 | 0.83 | 0.67 | 0.47 | 0.36 |
| **transport, copula by min. disagreement** | 0.42 | 0.38 | **0.17** | **0.15** | **0.16** |
| GP on $Y$ alone, 90% threshold | 0.29 | 0.24 | 0.23 | 0.20 | 0.21 |

For $\mathrm{CEVI} \geq 0.4$ the criterion recovers essentially all of the gap:
it attains the known-copula bound and beats a direct marginal fit. For
$\mathrm{CEVI} \leq 0.2$ it does not — it overshoots (median ratio 1.29 at
CEVI 0.1) and is more variable than simply fitting the marginal. That is the
ill-conditioned regime, and the boundary is now measured rather than asserted.

Two honest caveats. The chosen copula's disagreement sometimes falls *below* the
true copula's, so the criterion can reconcile the estimated conditionals better
than the truth does — visible as a consistent few-percent overshoot. And the
candidate set is Gaussian and the two Clayton rotations only; no extreme-value
copula (Gumbel, Hüsler–Reiss) was offered, and survival Clayton was chosen in
every single cell, which suggests the corner shape matters more than the family
label.

## 6. What the decomposition buys that a marginal fit cannot

1. **The assumption goes where it holds.** The conditional is exactly
   generalized Pareto; the marginal is not, anywhere reachable.

2. **A built-in consistency check, and it has teeth.** Disagreement between
   covariate values is a detectable internal contradiction with no analogue in a
   marginal fit, it correlates 0.53–0.74 with the actual error, and minimising
   it is enough to recover the known-copula bound where the problem is well
   conditioned. This is now the strongest claim in the file.

3. **It answers the question actually being asked.** For compound extremes: not
   only what the 200-year runoff is, but which combination of rainfall and
   snowmelt produces it.

## 7. Open questions

- The candidate set needs extreme-value copulas. Survival Clayton won every
  cell, which is suspicious — is the criterion selecting a corner or a family?
- Can the overshoot be removed? The criterion is minimised at a copula that
  over-reconciles noisy conditionals; some penalty or held-out version may fix it.
- Below $\mathrm{CEVI} = 0.4$ nothing works. Is that a fundamental limit, or an
  artefact of estimating the conditional shape freely when the conditional tail
  is nearly exponential?
- Two predictors means no scalar $F_X(x)$; the rank of the forest's predictive
  summary stands in. Is there a better one-dimensional reduction — and does the
  agreement criterion detect a bad one?

---

## References

- Coia, V., Joe, H. and Nolde, N. (2024). Copula-based conditional tail indices.
  *Journal of Multivariate Analysis* **201**, 105268.
  [doi:10.1016/j.jmva.2023.105268](https://doi.org/10.1016/j.jmva.2023.105268).
  The conditional tail index factorises into a copula-based conditional extreme
  value index and the marginal tail index; introduces a parametric family with
  a *non-constant* CEVI, and another that reduces a heavy-tailed response to a
  light tail on conditioning. This is the bookkeeping every claim above rests
  on, and the reason the copula is fitted along the upper edge rather than in
  the upper-right corner.
- Coia, V. (2017). *Forecasting of Nonlinear Extreme Quantiles Using Copula
  Models*. PhD dissertation, University of British Columbia.
  [cIRcle 1.0342941](https://open.library.ubc.ca/collections/24/items/1.0342941).
- `igcop`: Computational Tools for the IG and IGL Copula Families. R package,
  CRAN; documentation at [igcop.netlify.app](https://igcop.netlify.app/). The
  integrated gamma family is the natural next candidate for the agreement
  criterion in section 5, whose current candidate set contains only families
  with a constant CEVI — which is very likely why survival Clayton won every
  cell.

---

### Where the code is

| | |
|---|---|
| transport, ensembles, combination rules | `R/copula_transport.R` |
| the derived-marginal process | `R/transport_lab.R` |
| pipeline integration (off by default) | `R/dl_transport_marginal.R`, `scripts/5-runoff_marginals.r` |
| head-to-head against a direct fit | `scripts/experiments/transport-vs-marginal-fit.R` |
| choosing the copula by agreement | `scripts/experiments/copula-by-agreement.R` |
| earlier experiments | `scripts/experiments/{median-of-marginals,copula-transport-sweep,copula-transport-by-x}.R` |
| interactive version of the experiments | `apps/8-copula-transport-lab` |
