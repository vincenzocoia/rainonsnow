# Analyses

Every analysis models **hourly runoff at POT peaks as a function of a few
predictors**. An analysis is one folder, and is defined by three choices in
its `analysis.yaml`:

| Choice | Field | Where the options live |
|---|---|---|
| Which predictors | `predictors:` | `feature_registry()` in `R/features.R` — add wrangled predictors (e.g. available snowpack water) here |
| Which distributional learning model | `model: type:` | `dl_model_registry()` in `R/dl_fit.R` — `llqr` (one predictor), `rqforest`; a copula-transport wrapper would register here |
| How to query it | `queries:` | `R/queries.R` — event probability, trigger thresholds, likeliest drivers |

Shared upstream data (download, hourly table, POT peaks) is built once by
`scripts/1-3` into `derived/`; nothing in an analysis changes it.

## Index

| Analysis | Predictors | Model | Question | Status |
|---|---|---|---|---|
| [`rain-llqr`](rain-llqr/) | rainfall | `llqr` | Rain needed for a T-year runoff peak, ignoring snowmelt | active |
| [`rain-snowmelt-rqforest`](rain-snowmelt-rqforest/) | rainfall + snowmelt | `rqforest` | Rain needed for a T-year runoff peak, given the snowmelt rate | active |

Keep this table current: one row per folder. Mark finished or abandoned
analyses `status: archived` in their yaml rather than deleting them.

## Running

```bash
Rscript scripts/run-analysis.R rain-llqr          # all stages
Rscript scripts/run-analysis.R rain-llqr 4 5      # just queries + figures
```

| Stage | Script | Writes (`out/` unless noted) |
|---|---|---|
| 1 | `scripts/analysis/1-training.R` | `training.rds` (predictors at peaks), `event_levels.rds` |
| 2 | `scripts/analysis/2-fit.R` | `models.rds`, `predictions.rds`, `diagnostics.rds` |
| 3 | `scripts/analysis/3-marginals.R` | `return_levels.rds` |
| 4 | `scripts/analysis/4-queries.R` | `drivers.rds`, `event_probability.rds`, `triggers.rds`, `likeliest.rds` |
| 5 | `scripts/analysis/5-figures.R` | `figs/*.png` |

`out/` is not committed (regenerate it); `figs/` is.

Browse every analysis in one place with the explorer app:

```r
shiny::runApp("apps/explorer")
```

## Conventions that make analyses comparable

- **The T-year event is defined without the DL model.** Stage 1 fits a GP to
  each cell's POT peaks and converts return periods to runoff levels, so every
  analysis asks about the same event. Stage 3 checks whether the DL model's
  own implied marginal agrees.
- **Event probability** is `P(runoff > z_T | predictors)`, using predictive
  distributions with a grafted GP tail (`tail:` block), so it is not forced to
  zero beyond the largest training value.
- **Trigger threshold**: the smallest value of the target predictor (rain)
  at which the event probability reaches 10 / 50 / 90 %, for each value of the
  other predictors. "Not reached" means not within the observed range.
- **Likeliest drivers**: `f(x | runoff > z_T) ∝ P(runoff > z_T | x) f(x)`, with
  `f(x)` a parametric marginal (+ copula for two predictors) fitted to the
  drivers at peaks.

## Starting a new analysis

1. Copy an existing folder, rename it, edit `analysis.yaml`.
2. If it needs a new predictor, add it to `feature_registry()`.
3. `Rscript scripts/run-analysis.R <name>`.
4. Add a row above, and write what you learn in the folder's `README.md`.
