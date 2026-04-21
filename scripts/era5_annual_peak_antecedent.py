"""
From an ERA5-Land hourly NetCDF (same export as ``download_data.py`` / ``download_era5_range.py``):

1. Spatially aggregate runoff to a single hourly series (area mean by default).
2. For each water year or calendar year, find the hour of **annual maximum** runoff.
3. Save that hour's runoff and the **antecedent** hourly series of precipitation,
   snowfall, liquid rain (precip − snowfall, floored at zero), and snowmelt in the
   window leading up to and including the peak hour.

Inputs are local files only — no API keys. Run after you have downloaded NetCDF.

Example:
  python scripts/era5_annual_peak_antecedent.py \\
      --nc data/era5_land_hourly_alps.nc \\
      --window-hours 168 \\
      --water-year \\
      --out-dir outputs/era5_peaks
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Annual max runoff timing + antecedent met from ERA5-Land NetCDF."
    )
    p.add_argument("--nc", type=Path, required=True, help="Path to ERA5-Land hourly NetCDF.")
    p.add_argument(
        "--window-hours",
        type=int,
        default=168,
        help="Number of hours ending at peak (inclusive). Default 168 = 7 days.",
    )
    p.add_argument(
        "--water-year",
        action="store_true",
        help="Use water year (Oct–Sep, label = year of September end) instead of calendar year.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory for CSV outputs (default: <repo>/outputs/era5_peaks).",
    )
    return p.parse_args()


SPATIAL_DIM_NAMES = frozenset(
    {"latitude", "longitude", "lat", "lon", "x", "y", "band"}
)


def _time_dim_name(da: xr.DataArray) -> str:
    for d in da.dims:
        if d not in SPATIAL_DIM_NAMES:
            return d
    raise ValueError("Could not infer time dimension from runoff_hourly.")


def _spatial_dims(da: xr.DataArray) -> list[str]:
    tdim = _time_dim_name(da)
    return [d for d in da.dims if d != tdim]


def _year_labels(time: pd.DatetimeIndex, water_year: bool) -> np.ndarray:
    if not water_year:
        return time.year.values
    # Oct–Sep water year: label = calendar year of the September that closes the WY
    return np.where(time.month >= 10, time.year + 1, time.year)


def main() -> None:
    args = parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    out_dir = args.out_dir or (repo_root / "outputs" / "era5_peaks")
    out_dir.mkdir(parents=True, exist_ok=True)

    ds = xr.open_dataset(args.nc)
    runoff = ds["runoff_hourly"]
    tdim = _time_dim_name(runoff)
    spatial = _spatial_dims(runoff)
    runoff_mean = runoff.mean(dim=spatial)

    times = pd.to_datetime(runoff_mean[tdim].values)
    s = pd.Series(runoff_mean.values, index=times, name="runoff_m")
    years = _year_labels(s.index, args.water_year)

    peak_rows = []
    antecedent_rows = []

    precip = ds["total_precipitation_hourly"].mean(dim=spatial)
    snowfall = ds["snowfall_hourly"].mean(dim=spatial)
    snowmelt = ds["snowmelt_hourly"].mean(dim=spatial)

    for y in sorted(np.unique(years)):
        mask = years == y
        if not np.any(mask):
            continue
        sub = s.loc[mask]
        if sub.empty:
            continue
        peak_time = sub.idxmax()
        peak_val = float(sub.max())
        peak_rows.append(
            {
                "year": int(y),
                "peak_time_utc": peak_time.isoformat(),
                "runoff_spatial_mean_m_s": peak_val,
            }
        )

        delta = pd.Timedelta(hours=args.window_hours - 1)
        start = peak_time - delta
        tslice = slice(np.datetime64(start), np.datetime64(peak_time))

        p_win = precip.sel({tdim: tslice})
        sf_win = snowfall.sel({tdim: tslice})
        sm_win = snowmelt.sel({tdim: tslice})
        t_win = pd.to_datetime(p_win[tdim].values)

        rain = np.maximum(p_win.values - sf_win.values, 0.0)
        for i, t in enumerate(t_win):
            antecedent_rows.append(
                {
                    "year": int(y),
                    "time_utc": pd.Timestamp(t).isoformat(),
                    "hours_before_peak": float(
                        (peak_time - pd.Timestamp(t)).total_seconds() / 3600.0
                    ),
                    "total_precipitation_hourly": float(p_win.values[i]),
                    "snowfall_hourly": float(sf_win.values[i]),
                    "rainfall_hourly_proxy": float(rain[i]),
                    "snowmelt_hourly": float(sm_win.values[i]),
                }
            )

    peaks_df = pd.DataFrame(peak_rows).sort_values("year")
    ante_df = pd.DataFrame(antecedent_rows)

    ytag = "wy" if args.water_year else "cy"
    peaks_path = out_dir / f"annual_max_runoff_{ytag}.csv"
    ante_path = out_dir / f"antecedent_hourly_{ytag}_w{args.window_hours}.csv"
    peaks_df.to_csv(peaks_path, index=False)
    ante_df.to_csv(ante_path, index=False)
    print(f"Wrote {peaks_path}")
    print(f"Wrote {ante_path}")


if __name__ == "__main__":
    main()
