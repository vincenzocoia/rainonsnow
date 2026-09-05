# Notes on `probaverse/distionary`, branch `expectiles`

The expectile machinery on that branch still works after the supports change —
it installs clean against current dependencies, and `eval_expectile()` agrees
with an independent closed-form implementation to ~1e-9 for GEV, and to ~5e-11
for the normal, across shapes from xi = -0.2 to +0.4 and in the Gumbel limit.
`partial_moment()` is correct everywhere I tested it.

Two things are worth acting on.

## 1. Silently wrong expectiles for heavy tails (xi >~ 0.5)

`eval_expectile(dst_gev(0, 1, 0.7), at = p)` returns wrong values, with no
error — just a stream of muffled "Integration routine ... Returning NaN"
warnings:

| p | `eval_expectile` | correct |
|---|---|---|
| 0.60 | 5.690197 | 3.785170 |
| 0.65 | 5.690197 | 4.386845 |
| 0.70 | 5.690197 | 5.128520 |
| 0.75 | 6.083606 | 6.083606 |
| 0.80 | 7.390054 | 7.390054 |
| 0.90 | 22.760788 | 12.768850 |
| 0.99 | 91.043150 | 65.595430 |

The returned values are exactly 2, 8 and 32 times the mean — they are
bracket-expansion endpoints, not roots. The correct values satisfy the
identification equation `p E[(X-e)^+] = (1-p) E[(e-X)^+]` to 1e-15; the returned
ones do not.

**Cause.** Two independent things combine.

*(a) `distionary_integrate` is erratic on a heavy tail.* `partial_moment()`
evaluates `integrate(survival, lower = x, upper = Inf)`. For xi = 0.7 the
survival function decays like `t^(-1/0.7)`, and R's `integrate` fails for
particular lower limits, essentially at random — neighbouring values are fine:

```
phi(5.1285  ) = 1.712569     correct 1.71256857
phi(5.12852 ) = NaN          correct 1.71256643      <- 2e-5 further along
phi(5.6     ) = 1.664352     correct 1.66435205
phi(5.690197) = NaN          correct 1.65562625
```

*(b) `bracket_decreasing` in `src/numerics.h` treats a NaN as "not bracketed
yet".* The loop advances on `if (g(hi) <= 0.0) return true;`, and `NaN <= 0.0`
is `false`, so a single failed integration makes the bracket step straight over
the root:

```
seed = mean = 2.845, step = 2.845   ->  bracket candidate [2.845, 5.690]
g(5.690) = NaN, so `NaN <= 0` is false  ->  expand to [5.690, 11.380]
                                            (the root, 5.1285, is now outside)
newton_decreasing then finds g < 0 across the whole interval, shrinks hi onto
lo, and returns lo = 5.690197.
```

**Suggested fixes.**

- In `bracket_decreasing` and `newton_decreasing`, treat a non-finite `g` as a
  failure and return NaN, rather than letting it drive the search. A wrong
  number that looks reasonable is much worse than an NA.
- Make the partial moment robust rather than patching the solver around it. Two
  options, both reliable for heavy tails:
  - integrate in probability space, `phi(x) = int_{F(x)}^{1} (Q(u) - x) du`,
    which is a bounded interval; or
  - give families a closed form where one exists. For the GEV,
    substituting `z = (1 + xi (t - mu)/sigma)^(-1/xi)` and integrating by parts
    gives, for `xi < 1`, `xi != 0`,
    ```
    phi(x) = (sigma / xi) [ gammainc_lower(1 - xi, z_x) - (1 - e^{-z_x}) z_x^{-xi} ]
    ```
    with `gammainc_lower(s, a) = pgamma(a, s) * gamma(s)`. In R that is one
    line, exact to machine precision, and about 90x faster than integrating
    numerically (`R/gev.R` here; validated in `scripts/98-validate.R`).

## 2. Speed

`eval_expectile()` costs about 480 ms for 40 levels on a GEV, because each
Newton step integrates the survival function numerically from R. That is fine
interactively and prohibitive inside an optimiser — this experiment needed
roughly 10^8 expectile evaluations. With the closed-form partial moment the same
call costs about 2 ms in R and 0.13 ms in C++.

If expectiles are going to be usable for *fitting* and not just inspection, the
partial moment is the thing to make cheap.

## Branch state

`expectiles` is one commit on top of a `development` that is now 67 commits
behind, so it will need a merge. Nothing in the expectile code conflicted with
the supports work as far as this experiment exercised it.
