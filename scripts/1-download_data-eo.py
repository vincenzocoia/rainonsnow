"""
Download earth observation data for the Alps bbox to NetCDF (Google Earth Engine).
"""
from pathlib import Path
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
    return {"earth_engine": merged_ee, "download": merged_dl}


_cfg = load_download_config()
_PROJECT_ID = str(_cfg["earth_engine"]["project_id"])
_COLLECTION_ID = str(_cfg["download"]["collection_id"])
FIRST_YEAR = int(_cfg["download"]["first_year"])
LAST_YEAR = int(_cfg["download"]["last_year"])
selected_variables = list(_cfg["download"]["variables"])

# %%
# Initialize the Earth Engine module (project_id from inputs/data_specifications.yaml)
# ee.Authenticate()
ee.Initialize(project=_PROJECT_ID)


def download_data(year: int):
    ic = ee.ImageCollection(_COLLECTION_ID).filterDate(f"{year}-01-01", f"{year + 1}-01-01")
    alps = ee.Geometry.Rectangle([5.0, 43.0, 16.5, 48.5])
    ds = xr.open_dataset(
        ic,
        engine="ee",
        geometry=alps,
        backend_kwargs={"fast_time_slicing": True},
    )
    return ds


data_dir = repo_root / "data" / "eo"
data_dir.mkdir(parents=True, exist_ok=True)

years = range(FIRST_YEAR, LAST_YEAR + 1)
n_years = len(years)

# %%
# Download data for each year
for i, year in enumerate(years, start=1):
    output_path = data_dir / f"era5_land_hourly_alps_{year}.nc"
    print(f"\n[{i}/{n_years}] {year} -> {output_path.name}", flush=True)
    t0 = time.perf_counter()

    ds = download_data(year)
    ds = ds[selected_variables]

    if output_path.exists():
        output_path.unlink()
    ds.to_netcdf(output_path, engine="netcdf4")

    elapsed = time.perf_counter() - t0
    mib = output_path.stat().st_size // 1_048_576
    print(f"  done in {elapsed / 60.0:.2f} min (~{mib} MiB)", flush=True)


# %%
