# ESA satellite data options for the rain-on-snow analysis

Assessment of which ESA / Copernicus Earth-observation products can realistically enter
this pipeline, written against the code as it stands (`scripts/1-*` … `scripts/7-*`,
`inputs/*.yaml`). Compiled September 2026; product status changes, so the
"verify before committing" items at the end matter.

## Bottom line

Yes — there are several, and snow wetness is the right instinct. But none of them can
replace ERA5-Land as the *driver* of this model, and it would damage the analysis to try.
The defensible uses are, in order of value-per-effort:

1. **Validate the ERA5-Land snowmelt and snow-state fields** that currently act as model
   predictors, using Sentinel-1 wet-snow and Sentinel-2 snow-cover products. This is the
   weakest link in the whole chain and the cheapest thing to fix.
2. **Add an EO-derived snowpack-state covariate** (wet-snow fraction, fractional snow
   cover) to `dl_rqforest` and test whether it improves the conditional runoff
   distribution on the 2016–present overlap period.
3. **Build an observed rain-on-snow event catalogue** for the Alps (2016–present) and
   check it against the events the model flags as extreme.
4. **Static per-cell ancillary layers** (glacier cover, permafrost, land cover) from the
   ESA CCI — only useful if the model is ever pooled across cells.

## What the analysis actually needs

From `inputs/data_specifications.yaml` and `inputs/distributional_learning.yaml`:

| Requirement | Current source | Notes |
|---|---|---|
| Response: `runoff_hourly` | ERA5-Land | hourly, 1950–2025, 0.1° |
| Predictor: `rainfall_hourly` | ERA5-Land (`total_precipitation_hourly − snowfall_hourly`) | |
| Predictor: `snowmelt_hourly` | ERA5-Land | a *model diagnostic*, not an observation |
| Grid | 0.1° (~8 × 11 km at 47°N), Alps bbox | |
| Sampling | continuous hourly, no gaps | POT with q99 + 72 h declustering (`inputs/pot_metadata.yaml`) |

The extreme-value machinery (`scripts/3`, `scripts/5`, the GP tail graft) depends on a
long, gap-free hourly record. That is the constraint everything below runs into.

## The hard constraint

No ESA product gives hourly, gap-free, multi-decadal fields at this grid:

- Sentinel-1 based products start **September 2016** and sample every **6 days** nominally
  (better in the Alps when ascending and descending passes and overlapping relative
  orbits are combined, but still 1–3 day sampling at best, at two fixed local overpass
  times, not hourly).
- Sentinel-2 based products start 2016 and are **cloud-limited** — which is exactly the
  wrong failure mode, since rain-on-snow events happen under cloud.
- The long ESA CCI records that do go back to 1979–2000 are 1–25 km and, for SWE,
  explicitly mask Alpine terrain.

So: ERA5-Land stays as the driver. Satellite data enters as **evaluation**, as an
**additional covariate over the overlap period**, or as an **independent event record**.
Any proposal to ESA should say this plainly — it is a defensible scientific position, not
a concession.

## Recommended products

### 1. Copernicus SAR Wet Snow (SWS) — the "snow wetness" product

The direct match for the idea. Sentinel-1 C-band SAR, wet-snow extent in high-mountain
areas, **60 m**, Alps included, **September 2016 – present**, delivered per Sentinel-2
tile footprint (110 × 110 km). Part of the Copernicus Land Monitoring Service
High Resolution Water, Snow & Ice portfolio (HR-WSI; successor to HR-S&I).

Why it fits here: an ERA5-Land grid cell is ~8 × 11 km, so a 60 m wet-snow mask gives a
genuinely informative **sub-grid wet-snow fraction** — i.e. *how much of the cell had a
ripe, liquid-water-bearing snowpack* at the last overpass before a POT peak hour. That is
the physical state that separates a rain event that runs off from one that refreezes in
the pack, and it is precisely what `snowmelt_hourly` alone cannot tell you.

Two uses:

- **Validation.** Contingency table per cell and per season: ERA5-Land
  `snowmelt_hourly > ε` vs. observed wet snow. Gives melt-onset timing bias, which
  directly qualifies predictor quality in `scripts/4`.
- **Covariate.** `wet_snow_fraction` plus `hours_since_observation` added to
  `dl_rqforest` `xnames`, fitted on 2016–present only.

Caveats: binary wet/dry, not a quantitative liquid-water content — quantitative snow
wetness from SAR is still research-grade (2025 papers report MAE around 0.6–0.8 % LWC,
no operational product). Steep terrain causes layover/shadow gaps, which matters for your
current bbox (see below). Detection can miss gradual melt onset.

### 2. Fractional Snow Cover (FSC) and Gap-filled FSC (GFSC)

Sentinel-2 at 20 m (FSC) and the gap-filled, **daily** 60 m product (GFSC), Europe,
2016–present, same HR-WSI portfolio. GFSC is the more useful of the two here precisely
because it is gap-filled and daily — usable as a continuous covariate rather than a
sporadic observation.

Use it to **screen and weight POT events**: a "rain-on-snow" extreme in the current
pipeline is only inferred from ERA5-Land's own snow fields. FSC/GFSC tells you whether
snow was actually on the ground in that cell. Some fraction of the POT peaks currently
attributed to a rain+melt combination will turn out to have had little snow present; that
is a real finding, not just a data-quality note.

### 3. Sentinel-1 snow depth (KU Leuven C-SNOW, and the ALPSNOW / Digital Twin Alps line)

Sentinel-1 VH/VV change-detection snow depth over the European Alps at 500 m–1 km,
2017–present, reported R ≈ 0.87 and MAE ≈ 0.17 m against 743 Alpine in-situ sites; a
newer machine-learning version reaches 100 m. Use it to evaluate ERA5-Land
`snow_depth` / `snow_depth_water_equivalent`, which set the *melt capacity* — the ceiling
on how much water a rain-on-snow event can release.

ESA's own Alpine-regional projects are the most quotable option for an ESA-funded
proposal: **ALPSNOW** (snow extent, albedo, grain size, depth, SWE, melt area and
wetness, Alps-wide, 4 years) and **EO4Alps-snow**, whose snow model now runs in the
**Digital Twin Alps** project delivering daily **SWE, snow depth and snowmelt at 250 m
across the Alps**. A 250 m daily snowmelt field is the closest thing available to a direct
EO-informed replacement for the ERA5-Land `snowmelt_hourly` predictor — at daily, not
hourly, resolution, and over a short record. Being ESA project outputs, they also answer
the funder question directly.

### 4. ESA CCI Snow — snow cover fraction (the only long record)

Daily SCF at **1 km from MODIS (2000–2022)** and **5 km from AVHRR (1979–2022)**, with
per-pixel uncertainty and a canopy correction (snow viewable on canopy vs. on ground).
Coarse, but the only ESA snow record that overlaps a meaningful share of your 1950–2025
period. Use it for **snow-cover duration and climatology per cell back to the 1980s**, as
an independent check that ERA5-Land's Alpine snow seasonality is not drifting over the
record — relevant because a POT model fitted over 75 years assumes the driver fields are
homogeneous in time.

## Products that look relevant but are not usable here

| Product | Why not |
|---|---|
| **ESA CCI / GlobSnow SWE** (passive microwave, 1979–2023) | Complex terrain is **explicitly masked** — a mountain mask is applied from sub-grid elevation variability, because the 12.5–25 km footprint is incompatible with Alpine terrain. The Alps are a hole in this dataset. |
| **ESA CCI Soil Moisture** (combined, 0.25°, 1978–2024) | Coarse relative to your grid, and retrieval is **impossible over snow-covered and frozen ground** — masked exactly when you need it. Marginal value for *autumn antecedent* wetness before the snow season only. |
| Satellite precipitation | ESA does not produce a precipitation ECV. There is no ESA hourly satellite precipitation record for the Alps; radar/gauge products from national services, or EUMETSAT H SAF, are the alternatives, and neither is ESA. Precipitation stays with ERA5-Land. |
| **ESA CCI River Discharge** | Large-river discharge, not 0.1° grid-cell runoff. Not a substitute for the response variable. |
| **CCI Glaciers / Permafrost / High Resolution Land Cover** | Fine products, but static per-cell attributes are useless in the current design, where `dl_rqforest` is fitted **independently per cell** (`scripts/4-distributional_learning.r`). They only become predictors if the model is pooled or regionalised. Glacier mask is still worth having as a *diagnostic*: in high-Alpine cells, part of the ERA5-Land runoff signal is ice melt, not snowmelt. |

## Implementation paths

The pipeline is already Google Earth Engine based (`scripts/1-download_data-eo.py` via
`xee`), which splits the options cleanly:

**Path A — Sentinel-1 wet snow computed in Earth Engine (smallest change).**
`COPERNICUS/S1_GRD` is in the Earth Engine catalogue, already terrain-corrected and in dB.
Wet snow follows the standard Nagler–Rott change detection: ratio the wet-season
backscatter against a dry-snow/snow-free reference and threshold at about −3 dB (−2 dB
also used), with a layover/shadow and steep-slope mask. Reduce to the 0.1° ERA5-Land grid
as a per-cell wet fraction, and emit the same tabular shape `scripts/2` produces. New
script `scripts/1b-download_data-s1.py` plus a merge in `scripts/2`; add an `s1` block to
`inputs/data_specifications.yaml` alongside `download`.

Pros: reuses the existing Earth Engine auth, bbox and scale config. Cons: you are
reimplementing a classifier rather than using the validated operational one.

**Path B — download the operational HR-WSI products (better provenance).**
SWS / FSC / GFSC are not in the Earth Engine catalogue; they come from the Copernicus
Data Space Ecosystem (OData API) or WEkEO (HDA API), free of charge, as per-tile rasters
that need mosaicking and reprojection onto the 0.1° grid. More plumbing, but the product
is the official Copernicus one — which is the stronger claim in an ESA report.

**Recommended:** Path B for SWS and GFSC (provenance matters to the funder), with Path A
as a fallback if the archive access turns out to be awkward.

**Model-side design.** Keep the 1950–2025 ERA5-only model as the primary result. Fit a
parallel EO-informed model on the 2016/2017–present overlap, with
`xnames: [rainfall_hourly, snowmelt_hourly, wet_snow_fraction, fsc]`, and compare using
the skill machinery that already exists (`dl_skill_scores`, `dl_build_diagnostics`,
`apps/4b-distributional-learning-diagnostics`). That is a clean experimental design and
it uses ESA data to answer a real question rather than decoratively.

Be honest about the sample-size cost: q99 hourly POT with 72 h declustering over ~9
seasons leaves on the order of tens of peak events per cell. Enough for a skill
comparison on a subset of cells; not enough for the GP tail fits in `scripts/5`.

## Specific caveats for the current setup

- **The current bbox is a 2×2 test window** (`[10.75, 46.75, 10.95, 46.95]`, Ortler / Val
  Venosta). That is very steep terrain — check the fraction of each cell that is usable
  S1 backscatter after layover/shadow masking before designing anything around SWS. It
  may be low enough to force a different demonstration cell.
- **Sampling density is non-stationary.** Sentinel-1B failed in December 2021, degrading
  2022–2024 to a single-satellite 12-day baseline; Sentinel-1C (December 2024) and
  Sentinel-1D (November 2025) restore the two-satellite constellation, with the final
  C+D configuration reached around July 2026. A wet-snow covariate therefore has
  *time-varying missingness* across the record. That has to be handled explicitly —
  carry `hours_since_observation` as a companion covariate rather than silently
  interpolating.
- **Overpass timing vs. hourly peaks.** S1 passes the Alps at roughly 05:30 and 17:00
  local. A peak runoff hour can be up to ~72 h from the nearest observation; the covariate
  describes the pre-event snowpack state, not the state at the peak. Frame it that way.
- **Archive continuity.** HR-S&I production ended (catalogue retired through 2025–January
  2026) and HR-WSI took over; the SWS methodology was revised in 2024–2025 and the
  pre-2025 archive is being reprocessed for homogeneity, with publication to CDSE running
  through 2026. Pull a single consistent reprocessed version and store it locally — do
  not mix pre- and post-revision data in one time series.

## What to tell ESA

The honest framing, which is also the strongest one: satellite EO cannot drive a 75-year
hourly extreme-value model, but it is the only independent check on the snowpack-state
variables that model depends on, and it supplies a physically meaningful predictor
(sub-grid snow wetness) that reanalysis cannot provide. The project uses Copernicus
Sentinel-1 and Sentinel-2 products (SWS, GFSC) to validate and extend the ERA5-Land
drivers, and quantifies what observed snow wetness adds to the conditional distribution
of extreme runoff.

## Verify before committing

- Exact CDSE availability window for SWS and GFSC over the Alps *today*, given the
  ongoing reprocessing — and whether the pre-2025 reprocessed series is published yet.
- Usable (non-layover, non-shadow) S1 pixel fraction for the candidate cells.
- Licence and citation terms for the KU Leuven C-SNOW snow depth and for ALPSNOW /
  Digital Twin Alps products; the Copernicus HR-WSI products themselves are free and open.
- Whether the ESA–Polimi grant names specific products or missions it expects to see.

## Sources

- [Copernicus CLMS — SAR Wet Snow (SWS), 60 m, 2016–present](https://land.copernicus.eu/api/en/products/snow/high-resolution-sar-wet-snow)
- [Copernicus CLMS — Fractional Snow Cover (FSC), 20 m](https://land.copernicus.eu/api/en/products/snow/fractional-snow-cover)
- [Copernicus CLMS — Gap-filled Fractional Snow Cover (GFSC), 60 m, daily](https://land.copernicus.eu/api/en/products/snow/high-resolution-gap-filled-fractional-snow-cover)
- [CLMS — European snow and ice (HR-S&I) production concluded](https://land.copernicus.eu/en/production-updates/european-snow-and-ice-data-production-has-concluded-the-catalogue-will-remain-accessible)
- [CLMS — HR snow and ice archive 2016–2025 available until January 2026](https://land.copernicus.eu/en/production-updates/pan-european-high-resolution-snow-and-ice-data-archive-for-2016-2025-available-until-january-2026)
- [WEkEO — New Copernicus near-real-time products for snow and ice monitoring](https://wekeo.copernicus.eu/use-cases/new-copernicus-near-real-time-products-for-snow-and-ice-monitoring)
- [ESA Snow CCI project and data](https://climate.esa.int/en/projects/snow/Snow_data/)
- [Snow_cci Product User Guide v5.1](https://climate.esa.int/media/documents/Snow_cci_D4.3_PUG_v5.1.pdf)
- [GlobSnow v3.0 SWE dataset — complex-terrain masking](https://www.nature.com/articles/s41597-021-00939-2)
- [Evaluation of the Snow CCI snow-covered area product in a mountain SWE reanalysis](https://tc.copernicus.org/articles/19/2017/2025/)
- [ESA CCI Soil Moisture — project and data](https://climate.esa.int/en/projects/soil-moisture/data/)
- [ESA CCI Soil Moisture Product User Guide v09.1](https://climate.esa.int/documents/2960/ESA_CCI_SM_RD_D4.2_v2_Product_Users_Guide_v09.1_i1.0.pdf)
- [ESA CCI project list](https://climate.esa.int/en/projects/)
- [Lievens et al. (2022) — Sentinel-1 snow depth retrieval over the European Alps](https://tc.copernicus.org/articles/16/159/2022/)
- [Dunmire et al. — machine-learning snow depth across the European Alps from Sentinel-1](https://www.sciencedirect.com/science/article/pii/S003442572400395X)
- [ESA ALPSNOW project](https://eo4society.esa.int/projects/alpsnow/)
- [ESA EO4Alps-snow project](https://eo4society.esa.int/projects/eo4alps-snow/)
- [Snow cover products from Digital Twin Alps for alpine water management](https://eo4society.esa.int/2024/11/20/snow-cover-products-from-digital-twin-alps-for-alpine-water-management/)
- [Exploring how Sentinel-1 wet-snow maps can inform snowpack models](https://tc.copernicus.org/articles/18/5753/2024/)
- [Monitoring wet snow over an Alpine region using Sentinel-1](https://www.mdpi.com/2072-4292/13/3/381)
- [Wet and dry snow detection using Sentinel-1 SAR in mountainous areas](https://www.mdpi.com/2072-4292/11/8/895)
- [Snow wetness retrieval for snowmelt monitoring using dual-polarization SAR](https://www.tandfonline.com/doi/full/10.1080/10095020.2025.2575799)
- [Earth Engine — Sentinel-1 algorithms and COPERNICUS/S1_GRD](https://developers.google.com/earth-engine/guides/sentinel1)
- [Sentinel-1D user data opening and future plans](https://sentinels.copernicus.eu/-/sentinel-1d-user-data-opening-and-future-plans)
- [Sentinel-1 orbital reconfiguration dates](https://dataspace.copernicus.eu/news/2026-5-28-sentinel-1-orbital-reconfiguration-dates)
