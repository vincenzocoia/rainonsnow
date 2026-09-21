# Runoff given rainfall (local linear quantile regression)

**Question.** How much hourly rainfall does it take to produce a T-year hourly
runoff peak, ignoring snowmelt?

**Run.** 2 × 2 central-Alps grid, ERA5-Land 1950–2025, 2,566 POT peaks.
`llqr`, span 0.4, tricube, GP tail grafted at the median. Full run ≈ 4 min.

## Findings (2026-09-14)

- **Cell 3** (focus): rainfall explains 93% of runoff variance at peaks. Rain
  for a 50% chance of the event: 5.2 mm/h (2-yr), 6.7 (10-yr), 8.0 (50-yr),
  8.4 (100-yr). The 10–90% band is narrow, about 0.5–1 mm/h.
- **Cell 1** behaves like cell 3 (R² 91%): 5.1 mm/h (2-yr), 6.7 (10-yr).
- **Cells 2 and 4**: rainfall explains 1–3%. The 2-yr trigger is "not
  reached" in cell 2, and the 10-yr trigger in cell 4, because runoff there
  isn't driven by the hour's rain. Candidates for an accumulated or
  available-water predictor.
- **Tails end.** GP tails grafted onto the step distributions have negative
  shape in cells 1 and 3, so the event probability is exactly zero below
  ~5 mm/h (`figs/12-predictive-exceedance.png`). The DL-implied marginal
  flattens past ~70 years in those cells and overshoots in cell 4
  (`figs/02-return-levels.png`).

## Figures

- `10-llqr-quantile-curves`, `11-llqr-local-fit-{1,2,3}`: how the model works
  (local window at the 25th, 75th and 97th percentile of rain).
- `12-predictive-exceedance`: conditional runoff exceedance curves at six
  rainfall values.
- `20-event-probability`: P(event | rain), with 50% crossings.
- `21-likeliest-rainfall`: f(rain | event) against rain at all peaks.

## Open threads

- Span sensitivity (0.2 / 0.4 / 0.6) on the triggers.
- Out-of-sample calibration: diagnostics are in-sample.
