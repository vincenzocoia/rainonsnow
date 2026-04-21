"""
Download earth observation data for the Alps bbox to NetCDF (Google Earth Engine).
"""
from pathlib import Path
import time

import ee
import xarray as xr
import xee


COLLECTION_ID = "ECMWF/ERA5_LAND/HOURLY"

# Approximate coverage of ECMWF/ERA5_LAND/HOURLY on Earth Engine — edit when the catalog grows.
FIRST_YEAR = 1950
LAST_YEAR = 2025


# %%
# Initialize the Earth Engine module
# ee.Authenticate()
ee.Initialize(project="alps-data-explorer")


def download_data(year: int):
    ic = ee.ImageCollection(COLLECTION_ID).filterDate(f"{year}-01-01", f"{year + 1}-01-01")
    alps = ee.Geometry.Rectangle([5.0, 43.0, 16.5, 48.5])
    ds = xr.open_dataset(
        ic,
        engine="ee",
        geometry=alps,
        backend_kwargs={"fast_time_slicing": True},
    )
    return ds


# %%
# Subset to the following variables:
selected_variables = [
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


repo_root = Path(__file__).resolve().parent.parent
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
