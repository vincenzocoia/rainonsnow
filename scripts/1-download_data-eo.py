"""
Download earth observation data for the Alps bbox to NetCDF (Google Earth Engine).
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

import ee
import xarray as xr
import xee
import yaml


repo_root = Path(__file__).resolve().parent.parent
SPEC_PATH = repo_root / "inputs" / "data_specifications.yaml"

_DEFAULT_VARS = [
    "temperature_2m",
    "total_precipitation_hourly",
    "snowfall_hourly",
    "snowmelt_hourly",
    "snow_depth",
    "snow_depth_water_equivalent",
    "runoff_hourly",
    "surface_runoff_hourly",
    "sub_surface_runoff_hourly",
]


def load_download_config() -> dict:
    """Read inputs/data_specifications.yaml; fall back to embedded defaults."""
    base_ee = {"project_id": "alps-data-explorer"}
    base_dl = {
        "collection_id": "ECMWF/ERA5_LAND/HOURLY",
        "first_year": 1950,
        "last_year": 2025,
        "bbox": [5.0, 43.0, 16.5, 48.5],
        "scale": 0.1,
        "max_workers": 4,
        "variables": list(_DEFAULT_VARS),
    }
    if not SPEC_PATH.is_file():
        print(
            f"Note: {SPEC_PATH} missing; using defaults (see inputs/data_specifications.yaml).",
            flush=True,
        )
        return {"earth_engine": base_ee, "download": base_dl}
    with SPEC_PATH.open(encoding="utf-8") as f:
        spec = yaml.safe_load(f) or {}
    ee_cfg = spec.get("earth_engine") or {}
    dl_cfg = spec.get("download") or {}
    merged_ee = {**base_ee, **ee_cfg}
    merged_dl = {**base_dl, **dl_cfg}
    if merged_dl.get("variables") is None or merged_dl["variables"] == []:
        merged_dl["variables"] = base_dl["variables"]
    if merged_dl.get("bbox") is None:
        merged_dl["bbox"] = base_dl["bbox"]
    if merged_dl.get("max_workers") is None:
        merged_dl["max_workers"] = base_dl["max_workers"]
    if merged_dl.get("scale") is None:
        merged_dl["scale"] = base_dl["scale"]
    return {"earth_engine": merged_ee, "download": merged_dl}


_cfg = load_download_config()
_PROJECT_ID = str(_cfg["earth_engine"]["project_id"])
_COLLECTION_ID = str(_cfg["download"]["collection_id"])
FIRST_YEAR = int(_cfg["download"]["first_year"])
LAST_YEAR = int(_cfg["download"]["last_year"])
selected_variables = list(_cfg["download"]["variables"])
BBOX = [float(v) for v in _cfg["download"]["bbox"]]  # [lon_min, lat_min, lon_max, lat_max]
MAX_WORKERS = max(1, int(_cfg["download"]["max_workers"]))
SCALE = float(_cfg["download"]["scale"])  # output pixel size in degrees (ERA5-Land native = 0.1)

# %%
# Initialize the Earth Engine module (project_id from inputs/data_specifications.yaml)
# ee.Authenticate()
ee.Initialize(project=_PROJECT_ID)


def download_data(year: int):
    ic = ee.ImageCollection(_COLLECTION_ID).filterDate(f"{year}-01-01", f"{year + 1}-01-01")
    alps = ee.Geometry.Rectangle(BBOX)
    # `scale` (degrees) is required: without it xee picks a coarse default that
    # collapses a small box to a single pixel. 0.1 = ERA5-Land native resolution.
    ds = xr.open_dataset(
        ic,
        engine="ee",
        geometry=alps,
        scale=SCALE,
        backend_kwargs={"fast_time_slicing": True},
    )
    return ds


data_dir = repo_root / "derived" / "eo"
data_dir.mkdir(parents=True, exist_ok=True)

years = list(range(FIRST_YEAR, LAST_YEAR + 1))
n_years = len(years)


# A complete year's file is ~0.5 MiB even for a 2x2 grid; a partial/aborted
# write leaves a tiny stub. Treat anything below this as invalid/incomplete.
MIN_VALID_BYTES = 1024
MAX_RETRIES = 2  # attempts beyond the first, for transient EE errors


def download_year(year: int) -> dict:
    """Download and write one year's NetCDF. Returns a small status dict.

    Each Earth Engine year is an independent, network-bound request (~2 min of
    latency regardless of grid size), so years are fetched concurrently below.
    Resumable: a year whose file already exists and looks complete is skipped,
    so re-running after an interruption only fetches what is missing.
    """
    output_path = data_dir / f"era5_land_hourly_alps_{year}.nc"
    if output_path.exists() and output_path.stat().st_size >= MIN_VALID_BYTES:
        return {"year": year, "ok": True, "elapsed": 0.0,
                "kib": output_path.stat().st_size // 1024, "skipped": True}

    t0 = time.perf_counter()
    last_err = None
    for attempt in range(1 + MAX_RETRIES):
        try:
            ds = download_data(year)
            ds = ds[selected_variables]
            if output_path.exists():
                output_path.unlink()
            ds.to_netcdf(output_path, engine="netcdf4")
            return {
                "year": year,
                "ok": True,
                "elapsed": time.perf_counter() - t0,
                "kib": output_path.stat().st_size // 1024,
                "skipped": False,
            }
        except Exception as exc:  # keep other years going; retry transient errors
            last_err = exc
            if attempt < MAX_RETRIES:
                time.sleep(5 * (attempt + 1))
    return {"year": year, "ok": False, "elapsed": time.perf_counter() - t0,
            "error": repr(last_err)}


# %%
# Download years concurrently (I/O-bound EE requests overlap). Worker count from
# download.max_workers in inputs/data_specifications.yaml.
workers = min(MAX_WORKERS, n_years)
print(
    f"Downloading {n_years} year(s) [{FIRST_YEAR}-{LAST_YEAR}] with {workers} worker(s) "
    f"-> {data_dir}",
    flush=True,
)
t_all = time.perf_counter()
results = []
with ThreadPoolExecutor(max_workers=workers) as ex:
    futures = {ex.submit(download_year, y): y for y in years}
    done = 0
    for fut in as_completed(futures):
        res = fut.result()
        done += 1
        if res["ok"] and res.get("skipped"):
            print(
                f"[{done}/{n_years}] {res['year']} skipped (already downloaded, "
                f"~{res['kib']} KiB)",
                flush=True,
            )
        elif res["ok"]:
            print(
                f"[{done}/{n_years}] {res['year']} done in {res['elapsed'] / 60.0:.2f} min "
                f"(~{res['kib']} KiB)",
                flush=True,
            )
        else:
            print(
                f"[{done}/{n_years}] {res['year']} FAILED after {res['elapsed'] / 60.0:.2f} min: "
                f"{res['error']}",
                flush=True,
            )
        results.append(res)

failed = [r["year"] for r in results if not r["ok"]]
print(
    f"\nAll years finished in {(time.perf_counter() - t_all) / 60.0:.2f} min "
    f"({n_years - len(failed)}/{n_years} ok).",
    flush=True,
)
if failed:
    print(f"Failed years: {sorted(failed)}", flush=True)

# %%
