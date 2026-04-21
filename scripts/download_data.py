from pathlib import Path

import ee
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import xee


# %%
# Initialize the Earth Engine module
# ee.Authenticate()
ee.Initialize(project="alps-data-explorer")


# %%
# Download some ERA5 data
def download_data():
    ic = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY").filterDate(
        #"1992-10-05", "1993-03-31"
        "1992-01-01", "1992-12-31"
    )
    alps = ee.Geometry.Rectangle([5.0, 43.0, 16.5, 48.5])
    ds = xr.open_dataset(
        ic,
        engine="ee",
        geometry=alps,
    )
    return ds


ds = download_data()


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
ds = ds[selected_variables]


# %%
# Export to netCDF
repo_root = Path(__file__).resolve().parent.parent
data_dir = repo_root / "data"
data_dir.mkdir(parents=True, exist_ok=True)
output_path = data_dir / "era5_land_hourly_alps.nc"

ds.to_netcdf(output_path, engine="netcdf4")


# %%
