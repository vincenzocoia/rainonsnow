# Runoff given rainfall and snowmelt (quantile regression forest)

**Question.** For a T-year runoff peak, how much rain is enough when a given
amount of snowmelt is already happening? Snowmelt stands in for "water
available in the snowpack" until that predictor exists.

**Run.** Same grid and peaks as `rain-llqr`. `rqforest` (500 trees, mtry 1,
nodesize 23), GP tail grafted at the median, 50 × 50 query grid.
Full run ≈ 5 min.

## Findings (2026-09-14)

- **Cell 3, 2-year event**: melt lowers the rain needed, roughly one for one.
  50% chance at 6.6 mm/h rain with 0.31 mm/h melt, and 5.4 mm/h rain with
  1.0 mm/h melt. That matches a linear fit, where rain + melt explain 99% of
  runoff at peaks in cells 1 and 3.
- **5-year and rarer events cannot be resolved by the forest.** Only 2 peaks
  in cell 3 have more than 6 mm/h of rain, and a forest repeats its last
  leaves beyond the data. The chance of a 10-year peak never exceeds 65%, so
  most 50% thresholds are "not reached" and 90% is never reached. The event
  surfaces are blocky for the same reason.
- Cells 2 and 4: hourly rain and melt together explain 1–3% of runoff.

## Figures

- `20-trigger-rain-given-snowmelt`: rain needed against snowmelt, for 2-, 5-
  and 10-yr events at 10 / 50 / 90% chance.
- `21-event-surface-T{002,010,050}`: P(event) over (rain, melt).
- `22-likeliest-rain-given-snowmelt{,-T002}`: f(rain | event, melt).

## Open threads

- A smoother two-predictor learner, or copula transport, to get past the
  forest's ceiling for rarer events.
- Replace snowmelt with an available-water predictor
  (`feature_registry()`).
