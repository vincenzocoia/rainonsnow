"""
Derive Sentinel-1 wet-snow fraction per ERA5-Land grid cell (Google Earth Engine).

This is the analysis-ready companion to 1b-download_data-hrwsi.py. Rather than
downloading 60 m rasters and regridding them locally, it runs the classic
Nagler-Rott wet-snow change detection inside Earth Engine and returns, for every
Sentinel-1 acquisition, the *fraction of each 0.1 deg ERA5-Land cell classified as
wet snow*. That is the snowpack-preconditioning covariate: how much of the cell
held a ripe, liquid-water-bearing pack when the satellite passed over.

Method (Nagler & Rott 2000; Nagler et al. 2016):
  wet snow where  sigma0_obs / sigma0_ref  <  -3 dB
The reference is built per relative orbit and pass direction, as the mean
backscatter over mid-winter dry-snow months, so that observation geometry is
matched. Pixels are masked for permanent water and for local incidence angles
outside a usable range, which removes most radar layover and shadow.

Output (one CSV per year, resumable):
    derived/eo/s1_wetsnow/s1_wetsnow_<year>.csv
    columns: cell_id, x, y, datetime, relative_orbit, orbit_pass,
             wet_fraction, valid_pixels, total_pixels

`cell_id`, `x` and `y` follow the same 0.1 deg grid as
1-download_data-eo.py, so the result joins onto the ERA5-Land table from
2-tablify_spatial_eo.r on (x, y) -- but note the cell_id numbering there is
assigned in script 2, so join on coordinates, not on cell_id.

Caveats worth carrying into the modelling:
  * Sentinel-1 samples at ~05:30 and ~17:00 local, every 6-12 days per track.
    The value describes the pre-event snowpack state, not the state at the
    hourly runoff peak. Carry the time gap as a companion covariate.
  * Revisit density is not stationary: S1B failed in Dec 2021 (12-day baseline
    2022-2024), S1C arrived Dec 2024 and S1D Nov 2025.
  * In very steep terrain the usable pixel fraction can be small. Check
    `valid_pixels` / `total_pixels` before trusting a cell.
  * The detection is binary wet/dry, not a liquid water content.

Usage:
    earthengine authenticate          # once per machine
    uv run python scripts/1c-download_data-s1_wetsnow.py
    uv run python scripts/1c-download_data-s1_wetsnow.py --years 2019 2020
"""

from __future__ import annotations

import argparse
import csv
import math
from datetime import datetime, timezone
from pathlib import Path

import ee
import yaml

repo_root = Path(__file__).resolve().parent.parent
SPEC_PATH = repo_root / "inputs" / "data_specifications.yaml"

_DEFAULTS = {
    "first_year": 2017,
    "last_year": 2025,
    "months": [10, 11, 12, 1, 2, 3, 4, 5, 6],
    "polarisation": "VH",
    "threshold_db": -3.0,
    "reference_months": [12, 1, 2],
    "reference_percentile": 50,
    "lia_range": [20.0, 70.0],
    "reduce_scale": 60,
    "cells_per_request": 200,
    "output_dir": "derived/eo/s1_wetsnow",
}

# Sentinel-1 is right-looking; these are the nominal ground look azimuths
# (degrees clockwise from north) for the two pass directions.
LOOK_AZIMUTH = {"ASCENDING": 77.0, "DESCENDING": 283.0}


def load_config() -> dict:
    cfg = dict(_DEFAULTS)
    project_id = "alps-data-explorer"
    bbox = [10.75, 46.75, 10.95, 46.95]
    scale = 0.1
    if SPEC_PATH.is_file():
        with SPEC_PATH.open(encoding="utf-8") as f:
            spec = yaml.safe_load(f) or {}
        project_id = (spec.get("earth_engine") or {}).get("project_id", project_id)
        download = spec.get("download") or {}
        bbox = download.get("bbox", bbox)
        scale = download.get("scale", scale)
        for key, value in (spec.get("s1_wetsnow") or {}).items():
            if value is not None:
                cfg[key] = value
    else:
        print(f"Note: {SPEC_PATH} missing; using built-in defaults.", flush=True)
    cfg["project_id"] = project_id
    cfg["bbox"] = bbox
    cfg["grid_scale"] = scale
    return cfg


def build_cells(bbox: list[float], scale: float) -> list[dict]:
    """The 0.1 deg cell centres and polygons covering the bbox.

    Mirrors the grid that Earth Engine produces for the ERA5-Land export: cell
    centres sit on multiples of `scale`, which is what the bbox comment in
    inputs/data_specifications.yaml means by aligning box edges to x.x5.
    """
    lon_min, lat_min, lon_max, lat_max = bbox
    cells = []

    def centres(low: float, high: float) -> list[float]:
        first = math.floor(low / scale + 0.5) * scale
        out = []
        value = first
        while value < high:
            if value > low:
                out.append(round(value, 6))
            value = round(value + scale, 6)
        return out

    cell_id = 0
    for lat in centres(lat_min, lat_max):
        for lon in centres(lon_min, lon_max):
            cell_id += 1
            half = scale / 2.0
            cells.append(
                {
                    "cell_id": cell_id,
                    "x": lat,  # script 2 stores latitude in x and longitude in y
                    "y": lon,
                    "bounds": [
                        round(lon - half, 6),
                        round(lat - half, 6),
                        round(lon + half, 6),
                        round(lat + half, 6),
                    ],
                }
            )
    return cells


def cells_from_era5(nc_dir: Path, scale: float) -> list[dict] | None:
    """Build the cell grid from an ERA5-Land NetCDF written by script 1.

    Preferred over recomputing the grid: it guarantees the cell centres are
    byte-identical to the ones the ERA5 table carries, so the wet-snow CSV joins
    onto it on (x, y) without a half-cell offset. Returns None when no file is
    available, in which case the caller falls back to build_cells().
    """
    files = sorted(nc_dir.glob("era5_land_hourly_alps_*.nc"))
    if not files:
        return None
    try:
        import xarray as xr
    except ImportError:
        return None
    try:
        ds = xr.open_dataset(files[0])
    except Exception as err:  # noqa: BLE001 - fall back rather than fail the run
        print(f"Note: could not read {files[0].name} ({err}); computing the grid instead.", flush=True)
        return None

    with ds:
        lon_name = next((n for n in ("lon", "longitude", "X", "x") if n in ds.coords), None)
        lat_name = next((n for n in ("lat", "latitude", "Y", "y") if n in ds.coords), None)
        if lon_name is None or lat_name is None:
            print(f"Note: no lon/lat coords in {files[0].name}; computing the grid instead.", flush=True)
            return None
        lons = [round(float(v), 6) for v in ds[lon_name].values]
        lats = [round(float(v), 6) for v in ds[lat_name].values]

    half = scale / 2.0
    cells = []
    cell_id = 0
    for lat in sorted(lats):
        for lon in sorted(lons):
            cell_id += 1
            cells.append(
                {
                    "cell_id": cell_id,
                    "x": lat,
                    "y": lon,
                    "bounds": [
                        round(lon - half, 6),
                        round(lat - half, 6),
                        round(lon + half, 6),
                        round(lat + half, 6),
                    ],
                }
            )
    print(f"Grid taken from {files[0].name}: {len(cells)} cells", flush=True)
    return cells


def cells_to_fc(cells: list[dict]) -> ee.FeatureCollection:
    features = [
        ee.Feature(
            ee.Geometry.Rectangle(cell["bounds"], "EPSG:4326", False),
            {"cell_id": cell["cell_id"], "x": cell["x"], "y": cell["y"]},
        )
        for cell in cells
    ]
    return ee.FeatureCollection(features)


def local_incidence_angle(image: ee.Image, orbit_pass: str) -> ee.Image:
    """Terrain-corrected incidence angle, from the DEM and the nominal look azimuth."""
    dem = ee.ImageCollection("COPERNICUS/DEM/GLO30").select("DEM").mosaic()
    terrain = ee.Algorithms.Terrain(dem)
    slope = terrain.select("slope").multiply(math.pi / 180.0)
    aspect = terrain.select("aspect").multiply(math.pi / 180.0)
    theta = image.select("angle").multiply(math.pi / 180.0)
    look = ee.Number(LOOK_AZIMUTH.get(orbit_pass, 77.0)).multiply(math.pi / 180.0)
    # cos(LIA) = cos(theta) cos(slope) + sin(theta) sin(slope) cos(look - aspect)
    cos_lia = (
        theta.cos().multiply(slope.cos())
        .add(theta.sin().multiply(slope.sin()).multiply(aspect.subtract(look).cos()))
    )
    return cos_lia.clamp(-1, 1).acos().multiply(180.0 / math.pi).rename("lia")


def base_collection(cfg: dict, region: ee.Geometry) -> ee.ImageCollection:
    pol = cfg["polarisation"]
    return (
        ee.ImageCollection("COPERNICUS/S1_GRD")
        .filterBounds(region)
        .filter(ee.Filter.eq("instrumentMode", "IW"))
        .filter(ee.Filter.listContains("transmitterReceiverPolarisation", pol))
        .select([pol])
    )


def reference_image(cfg: dict, region: ee.Geometry, relative_orbit: int, orbit_pass: str) -> ee.Image:
    """Dry-snow reference backscatter (linear power) for one track."""
    pol = cfg["polarisation"]
    # An Or of single-month filters, so a winter list that wraps the year end
    # (e.g. [12, 1, 2]) works the same as a contiguous one.
    months = ee.Filter.Or(
        *[ee.Filter.calendarRange(m, m, "month") for m in cfg["reference_months"]]
    )
    reference = (
        base_collection(cfg, region)
        .filter(ee.Filter.eq("relativeOrbitNumber_start", relative_orbit))
        .filter(ee.Filter.eq("orbitProperties_pass", orbit_pass))
        .filter(months)
    )
    # Average in linear power, not in dB: the -3 dB test is a power ratio.
    linear = reference.map(lambda img: ee.Image(10).pow(img.divide(10.0)).rename(pol))
    return linear.reduce(ee.Reducer.percentile([cfg["reference_percentile"]])).rename("ref")


def wet_snow_image(cfg: dict, image: ee.Image, reference: ee.Image, orbit_pass: str) -> ee.Image:
    """Binary wet-snow mask plus a validity mask, for one acquisition."""
    pol = cfg["polarisation"]
    observed = ee.Image(10).pow(image.select(pol).divide(10.0))
    ratio_db = observed.divide(reference).log10().multiply(10.0)

    lia = local_incidence_angle(image, orbit_pass)
    lia_ok = lia.gte(cfg["lia_range"][0]).And(lia.lte(cfg["lia_range"][1]))

    water = (
        ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence").unmask(0).lt(50)
    )

    valid = image.select(pol).mask().And(reference.mask()).And(lia_ok).And(water)
    wet = ratio_db.lt(cfg["threshold_db"]).rename("wet").updateMask(valid)
    return wet


def pixel_totals(cfg: dict, cells_fc: ee.FeatureCollection) -> dict:
    """Unmasked pixel count per cell at the reduction scale.

    The denominator for `valid_pixels`: it depends only on the cells and the
    scale, so it is computed once per cell chunk rather than once per scene.
    """
    counts = ee.Image(1).rename("wet").reduceRegions(
        collection=cells_fc,
        reducer=ee.Reducer.count(),
        scale=cfg["reduce_scale"],
    )
    return {
        f["properties"]["cell_id"]: f["properties"].get("count", 0)
        for f in counts.getInfo()["features"]
    }


def rows_for_month(
    cfg: dict,
    region: ee.Geometry,
    cells_fc: ee.FeatureCollection,
    total_by_cell: dict,
    year: int,
    month: int,
) -> list[dict]:
    """Per-cell wet fraction for every acquisition in one month."""
    start = ee.Date.fromYMD(year, month, 1)
    end = start.advance(1, "month")
    collection = base_collection(cfg, region).filterDate(start, end)

    try:
        scenes = collection.toList(collection.size()).getInfo()
    except ee.EEException as err:
        print(f"    ! {year}-{month:02d}: could not list scenes ({err})", flush=True)
        return []
    if not scenes:
        return []

    rows: list[dict] = []
    references: dict[tuple[int, str], ee.Image] = {}

    for scene in scenes:
        props = scene.get("properties", {})
        relative_orbit = props.get("relativeOrbitNumber_start")
        orbit_pass = props.get("orbitProperties_pass")
        if relative_orbit is None or orbit_pass is None:
            continue
        key = (int(relative_orbit), str(orbit_pass))
        if key not in references:
            references[key] = reference_image(cfg, region, key[0], key[1])

        image = ee.Image(scene["id"])
        wet = wet_snow_image(cfg, image, references[key], orbit_pass)

        stats = wet.reduceRegions(
            collection=cells_fc,
            reducer=ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True),
            scale=cfg["reduce_scale"],
        )
        try:
            stats_info = stats.getInfo()["features"]
        except ee.EEException as err:
            print(f"    ! {scene['id']}: {err}", flush=True)
            continue

        timestamp = props.get("system:time_start")
        stamp = (
            datetime.fromtimestamp(timestamp / 1000.0, tz=timezone.utc).strftime(
                "%Y-%m-%dT%H:%M:%S"
            )
            if timestamp is not None
            else ""
        )
        for feature in stats_info:
            prop = feature["properties"]
            rows.append(
                {
                    "cell_id": prop["cell_id"],
                    "x": prop["x"],
                    "y": prop["y"],
                    "datetime": stamp,
                    "relative_orbit": key[0],
                    "orbit_pass": key[1],
                    "wet_fraction": prop.get("mean"),
                    "valid_pixels": prop.get("count", 0),
                    "total_pixels": total_by_cell.get(prop["cell_id"], 0),
                }
            )
    return rows


def write_year(path: Path, rows: list[dict]) -> None:
    fields = [
        "cell_id",
        "x",
        "y",
        "datetime",
        "relative_orbit",
        "orbit_pass",
        "wet_fraction",
        "valid_pixels",
        "total_pixels",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    cfg = load_config()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--years", nargs="+", type=int, help="override the year range")
    parser.add_argument("--overwrite", action="store_true", help="redo years already on disk")
    parser.add_argument(
        "--compute-grid",
        action="store_true",
        help="derive cells from the bbox instead of reading derived/eo/*.nc",
    )
    args = parser.parse_args()

    ee.Initialize(project=cfg["project_id"])

    bbox = cfg["bbox"]
    region = ee.Geometry.Rectangle(bbox, "EPSG:4326", False)
    cells = None
    if not args.compute_grid:
        cells = cells_from_era5(repo_root / "derived" / "eo", cfg["grid_scale"])
    if cells is None:
        cells = build_cells(bbox, cfg["grid_scale"])
        print(
            f"{len(cells)} grid cells computed from bbox {bbox}. No ERA5-Land NetCDF "
            "was read, so check that these centres match the ones in "
            "derived/era5_land_hourly_alps_all.rds before joining on (x, y).",
            flush=True,
        )
    if not cells:
        raise SystemExit(f"No {cfg['grid_scale']} deg cell centres fall inside bbox {bbox}.")

    years = args.years or list(range(int(cfg["first_year"]), int(cfg["last_year"]) + 1))
    months = sorted(set(cfg["months"])) if cfg["months"] else list(range(1, 13))
    out_dir = repo_root / cfg["output_dir"]

    chunk = int(cfg["cells_per_request"])
    if len(cells) > chunk:
        print(
            f"Note: {len(cells)} cells exceeds cells_per_request={chunk}; requests are "
            "chunked. For the full Alps grid this script will be slow -- consider "
            "ee.batch.Export.table.toDrive instead.",
            flush=True,
        )

    for year in years:
        out_path = out_dir / f"s1_wetsnow_{year}.csv"
        if out_path.is_file() and not args.overwrite:
            print(f"{year}: already present, skipping ({out_path.name})", flush=True)
            continue
        print(f"{year}: processing", flush=True)
        rows: list[dict] = []
        totals_cache: dict[tuple[int, int], dict] = {}
        for month in months:
            for start in range(0, len(cells), chunk):
                subset = cells[start : start + chunk]
                subset_fc = cells_to_fc(subset)
                key = (start, len(subset))
                if key not in totals_cache:
                    totals_cache[key] = pixel_totals(cfg, subset_fc)
                rows.extend(
                    rows_for_month(
                        cfg, region, subset_fc, totals_cache[key], year, month
                    )
                )
            print(f"  {year}-{month:02d}: {len(rows)} rows so far", flush=True)
        if rows:
            write_year(out_path, rows)
            print(f"{year}: wrote {len(rows)} rows to {out_path}", flush=True)
        else:
            print(f"{year}: no acquisitions matched, nothing written", flush=True)


if __name__ == "__main__":
    main()
