# Superseded

Code kept for reference after the `analyses/` restructure replaced it. Nothing
here is maintained, wired into the current pipeline, or expected to run against
the repository as it stands today. It is here to be *read*; git has the history,
but that is harder to browse than a directory.

## `apps/` — the pre-restructure Shiny apps

Eight explorers written against the old numbered pipeline (`scripts/4`–`7`,
`inputs/distributional_learning.yaml`, `inputs/rain_snow_joint_model.yaml`),
all of which the restructure removed. They read `derived/*.rds` files that the
current `Rscript scripts/run-analysis.R` no longer produces under those names —
per-analysis outputs now live in `analyses/<name>/out/`.

**To see them as they were, check out the tag:**

```bash
git checkout apps-0.1.0
```

That tag is the last commit where these apps and the pipeline feeding them
existed together (`7656717`, 2026-09-07, on the now-merged
`claude/rain-on-snow-stats-vc7ouw`). The app files there are byte-identical to
the copies in this folder — the tag adds the scripts, configs and package state
they expect. Running them still needs the `derived/` artifacts, which are
gitignored and must be regenerated.

| App | What it showed | Needed |
|---|---|---|
| `4a-distributional-learning-fit` | Fit and tune the distributional-learning model, one cell at a time; wrote `inputs/distributional_learning.yaml` | POT peaks |
| `4b-distributional-learning-diagnostics` | Full-grid diagnostics — P-P calibration, median quantile skill across τ | POT peaks |
| `5-runoff-marginals-explorer` | Mixture of hourly predictive distributions on a POT event axis | Peaks, DL return levels |
| `5b-tail-shape-explorer` | Per-cell GPD tail shape, each cell fitted on its own data, click to select | Mixture tails, tail summary |
| `6-joint-rain-snow-explorer` | Joint rainfall–snowmelt structure per cell | Hourly table, joint model |
| `7-rain-snow-given-runoff` | Conditional rain–snow structure given extreme runoff | Peaks, rqforest models, joint model, mixture tails |
| `8-copula-transport-lab` | Recovering a marginal tail from conditional tails, scored against a known truth | **nothing** — fully simulated |

`dl_shared.R` and `ros_theme.R` are the helpers they sourced.

`8-copula-transport-lab` is the exception worth knowing about: it reads no
`derived/` data at all — everything on screen is simulated from a
data-generating process chosen in the sidebar, so the truth is known exactly.
It is the interactive companion to `notes/marginal-from-conditionals.md`, and of
everything here it is the most likely to still run.

## What replaced them

`apps/explorer` — a single hub app driven by the same `analyses/<name>/` outputs
as the rest of the pipeline, sharing plot helpers with `scripts/lib/plots.R`.
`apps/3-pot-explorer` was not superseded and is still live.

## Caveats

These were research tools, not products. They had rough edges in their prime —
the last three commits to touch them were bug fixes for panels that errored
outright and axes reading longitude as latitude — and there is no claim that a
given panel worked on a given day. **Their functionality has not been
re-tested**, and re-testing it is not planned. Treat them as a record of what
was explored and how it was framed, not as working software.
