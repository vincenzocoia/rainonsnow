"""
Download ERA5-Land hourly fields for the Alps bbox to NetCDF (Google Earth Engine).

By default this uses the **full temporal coverage** of ``ECMWF/ERA5_LAND/HOURLY`` on
Earth Engine (queried at runtime) and writes **one NetCDF per calendar decade**
(chunks are half-open [chunk_start, chunk_end) aligned to Jan 1 of years 1960, 1970, …).

Optional ``--start`` and ``--end`` restrict the span (both required together); output
is still split by decade within that span.

Prerequisites (manual / one-time):
  - Google Cloud project with Earth Engine enabled
  - ``ee.Authenticate()`` once on this machine
  - ``ee.Initialize(project=...)`` with your project id

No Copernicus CDS key is required for this path; GEE uses its own auth.

Examples:
  python scripts/download_era5_range.py --project alps-data-explorer

  python scripts/download_era5_range.py --project alps-data-explorer \\
      --out-dir data --name-prefix era5_land_hourly_alps

  python scripts/download_era5_range.py --start 1990-01-01 --end 2000-01-01 \\
      --project alps-data-explorer

Runtime (very rough): each decade for this bbox and variable set is often on the
order of **tens of minutes to a few hours** depending on Earth Engine load, your
network, and whether data are warm in the cache; the full catalog (~7 decades)
can total **many hours** and may occasionally hit memory or quota limits—in that
case, narrow the date range or pull fewer variables in a forked script.

Within each decade file, data are written **year-by-year** to small segment
NetCDFs, then **merged** into one file (works on older xarray that lack
``Dataset.to_netcdf(..., append_dim=...)``). Each year stays under xee’s ~10k
hourly image soft limit. **fast_time_slicing** is enabled (see xee docs).

While each year is written, a **heartbeat** line prints every 60 seconds so a
long single-year export still shows the process is alive.
"""

from __future__ import annotations

import argparse
import logging
import threading
import time
from pathlib import Path

import ee
import pandas as pd
import xarray as xr
import xee

COLLECTION_ID = "ECMWF/ERA5_LAND/HOURLY"

def _configure_logging() -> None:
    """Drop repetitive xee warnings when the hourly series exceeds 10k steps."""
    class _SuppressXee10kWarning(logging.Filter):
        def filter(self, record: logging.LogRecord) -> bool:
            return "beyond 10000 images" not in record.getMessage()

    logging.getLogger().addFilter(_SuppressXee10kWarning())


SELECTED_VARIABLES = [
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="ERA5-Land hourly download via Earth Engine + xee (decade NetCDFs)."
    )
    p.add_argument(
        "--start",
        default=None,
        help="Optional start date (YYYY-MM-DD). If set, --end is required.",
    )
    p.add_argument(
        "--end",
        default=None,
        help="Optional end date (YYYY-MM-DD), passed to filterDate as exclusive upper bound.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Directory for NetCDF files (default: <repo>/data).",
    )
    p.add_argument(
        "--name-prefix",
        default="era5_land_hourly_alps",
        help="Filename prefix before the date span (default: era5_land_hourly_alps).",
    )
    p.add_argument(
        "--project",
        default="alps-data-explorer",
        help="Google Cloud / Earth Engine project id for ee.Initialize.",
    )
    p.add_argument(
        "--west", type=float, default=5.0, help="BBox west longitude (deg E)."
    )
    p.add_argument(
        "--south", type=float, default=43.0, help="BBox south latitude (deg N)."
    )
    p.add_argument(
        "--east", type=float, default=16.5, help="BBox east longitude (deg E)."
    )
    p.add_argument(
        "--north", type=float, default=48.5, help="BBox north latitude (deg N)."
    )
    return p.parse_args()


def _validate_optional_range(start: str | None, end: str | None) -> None:
    if (start is None) ^ (end is None):
        raise SystemExit("Provide both --start and --end, or neither for full catalog span.")


def collection_time_range(collection_id: str) -> tuple[pd.Timestamp, pd.Timestamp]:
    """Return (t_start, t_end) in UTC; chunking uses half-open intervals up to t_end."""
    ic = ee.ImageCollection(collection_id)
    t0_ms = ic.aggregate_min("system:time_start").getInfo()
    t1_ms = ic.aggregate_max("system:time_end").getInfo()
    if t1_ms is None:
        last = ic.sort("system:time_start", False).limit(1).first()
        t1_ms = last.get("system:time_end").getInfo()
        if t1_ms is None:
            t0_last = last.get("system:time_start").getInfo()
            t1_ms = int(t0_last) + 3_600_000
    t0 = pd.to_datetime(int(t0_ms), unit="ms", utc=True)
    t1 = pd.to_datetime(int(t1_ms), unit="ms", utc=True)
    return t0, t1


def _utc_ts(ts: pd.Timestamp) -> pd.Timestamp:
    if ts.tzinfo is None:
        return ts.tz_localize("UTC")
    return ts.tz_convert("UTC")


def _gee_date_str(ts: pd.Timestamp) -> str:
    ts = _utc_ts(ts)
    return ts.strftime("%Y-%m-%dT%H:%M:%SZ")


def _file_tag(ts: pd.Timestamp) -> str:
    return _utc_ts(ts).strftime("%Y%m%dT%H%M%SZ")


def year_intervals(a: pd.Timestamp, b: pd.Timestamp) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Half-open sub-intervals [ya, yb) covering [a, b) aligned to calendar-year splits."""
    a = _utc_ts(a)
    b = _utc_ts(b)
    spans: list[tuple[pd.Timestamp, pd.Timestamp]] = []
    cur = a
    while cur < b:
        nxt = min(pd.Timestamp(year=cur.year + 1, month=1, day=1, tz="UTC"), b)
        if nxt <= cur:
            break
        spans.append((cur, nxt))
        cur = nxt
    return spans


def _infer_time_dimension(ds: xr.Dataset) -> str:
    for name in ("time", "valid_time"):
        if name in ds.dims:
            return name
    spatial = frozenset({"latitude", "longitude", "lat", "lon", "x", "y"})
    for d in ds.dims:
        if d not in spatial:
            return d
    raise ValueError("Could not infer time dimension for NetCDF append.")


def _heartbeat_thread(interval_s: float) -> tuple[threading.Event, threading.Thread]:
    stop = threading.Event()

    def run() -> None:
        while not stop.wait(interval_s):
            print(f"    … still writing ({time.strftime('%H:%M:%S')})", flush=True)

    th = threading.Thread(target=run, name="era5-download-heartbeat", daemon=True)
    th.start()
    return stop, th


def decade_chunks(t0: pd.Timestamp, t1: pd.Timestamp) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """
    Half-open chunks [a, b) with b the earlier of the next Jan 1 on a multiple-of-10
    year boundary after ``a``, and ``t1`` (exclusive catalog end).
    """
    t0 = _utc_ts(t0)
    t1 = _utc_ts(t1)

    chunks: list[tuple[pd.Timestamp, pd.Timestamp]] = []
    cur = t0
    while cur < t1:
        decade_end_year = ((cur.year // 10) + 1) * 10
        boundary = pd.Timestamp(year=decade_end_year, month=1, day=1, tz="UTC")
        nxt = min(boundary, t1)
        if nxt <= cur:
            break
        chunks.append((cur, nxt))
        cur = nxt
    return chunks


def main() -> None:
    args = parse_args()
    _validate_optional_range(args.start, args.end)
    _configure_logging()

    repo_root = Path(__file__).resolve().parent.parent
    out_dir = (args.out_dir or (repo_root / "data")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    ee.Initialize(project=args.project)

    if args.start is None:
        t0, t1 = collection_time_range(COLLECTION_ID)
        print(
            f"Catalog {COLLECTION_ID} span (UTC): {t0.isoformat()} .. {t1.isoformat()} "
            "(end time is last image end; chunks are half-open [a, b))"
        )
    else:
        t0 = pd.Timestamp(args.start, tz="UTC")
        t1 = pd.Timestamp(args.end, tz="UTC")

    chunks = decade_chunks(t0, t1)
    if not chunks:
        raise SystemExit("No time chunks to download (empty range?).")

    alps = ee.Geometry.Rectangle([args.west, args.south, args.east, args.north])

    for a, b in chunks:
        fa, fb = _gee_date_str(a), _gee_date_str(b)
        out = out_dir / f"{args.name_prefix}_{_file_tag(a)}_{_file_tag(b)}.nc"
        print(f"Decade file {fa} .. {fb} -> {out.name}", flush=True)
        t_chunk0 = time.perf_counter()
        if out.exists():
            print(f"  Replacing existing file ({out.stat().st_size // 1_048_576} MiB)", flush=True)
            out.unlink()

        for stale in sorted(out_dir.glob(f"{out.stem}__seg*.nc")):
            stale.unlink()

        years = year_intervals(a, b)
        n_years = len(years)
        part_paths: list[Path] = []

        for i, (ya, yb) in enumerate(years, start=1):
            yfa, yfb = _gee_date_str(ya), _gee_date_str(yb)
            part = out_dir / f"{out.stem}__seg{i:02d}.nc"
            part_paths.append(part)
            print(f"  [{i}/{n_years}] year block {yfa} .. {yfb} -> {part.name}", flush=True)
            t_y0 = time.perf_counter()
            ic = ee.ImageCollection(COLLECTION_ID).filterDate(yfa, yfb)
            ds = xr.open_dataset(
                ic,
                engine="ee",
                geometry=alps,
                backend_kwargs={"fast_time_slicing": True},
            )
            ds = ds[SELECTED_VARIABLES]
            tdim = _infer_time_dimension(ds)
            nt = ds.sizes.get(tdim)
            if isinstance(nt, int) and nt > 0:
                print(f"      opening: {nt} hourly steps along {tdim!r}", flush=True)

            hb_stop, hb_thread = _heartbeat_thread(60.0)
            try:
                if part.exists():
                    part.unlink()
                ds.to_netcdf(part, engine="netcdf4")
            finally:
                hb_stop.set()
                hb_thread.join(timeout=5.0)

            dt_y = time.perf_counter() - t_y0
            sz_mib = part.stat().st_size // 1_048_576
            print(f"      segment written in {dt_y / 60.0:.2f} min (~{sz_mib} MiB)", flush=True)

        print(f"  Merging {len(part_paths)} segments -> {out.name}", flush=True)
        hb_stop, hb_thread = _heartbeat_thread(60.0)
        merge_ok = False
        try:
            paths = [str(p) for p in part_paths]
            try:
                merged = xr.open_mfdataset(
                    paths,
                    combine="by_coords",
                    engine="netcdf4",
                )
            except TypeError:
                merged = xr.open_mfdataset(paths, concat_dim="time", engine="netcdf4")
            try:
                merged.to_netcdf(out, engine="netcdf4")
                merge_ok = True
            finally:
                merged.close()
        finally:
            hb_stop.set()
            hb_thread.join(timeout=5.0)

        if merge_ok:
            for p in part_paths:
                if p.exists():
                    p.unlink()

        elapsed = time.perf_counter() - t_chunk0
        sz_final = out.stat().st_size // 1_048_576
        print(
            f"  Finished decade {out.name} (~{sz_final} MiB) in {elapsed / 60.0:.1f} min total",
            flush=True,
        )


if __name__ == "__main__":
    main()
