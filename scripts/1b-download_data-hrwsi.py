"""
Download Copernicus HR-WSI snow products (wet snow, snow cover) for the Alps bbox.

Companion to 1-download_data-eo.py (ERA5-Land). Where ERA5-Land supplies the
hourly driver fields, these Sentinel-1/Sentinel-2 products supply the observed
*snowpack preconditioning* state -- whether the pack was wet and how much of the
cell was snow covered -- at satellite revisit frequency from September 2016.

Products (see inputs/data_specifications.yaml -> hrwsi.product_types):
  SWS   SAR Wet Snow, 60 m, Sentinel-1, high-mountain areas. The wet-snow layer.
  GFSC  Gap-filled Fractional Snow Cover, 60 m, daily, S1+S2. Snow-covered fraction.
  WDS   Wet/Dry Snow, 20 m, Sentinel-2 (cloud limited).
  FSC   Fractional Snow Cover, 20 m, Sentinel-2 (cloud limited).

Access is anonymous: the Copernicus Land Monitoring Service publishes HR-WSI on a
CloudFerro S3 bucket with public read keys. The endpoint, bucket and keys below
are the ones shipped in the official EEA client
(https://github.com/eea/clms-hrwsi-api-client-python) and can be overridden via
the HRWSI_S3_* environment variables or the YAML config.

Usage:
    uv run python scripts/1b-download_data-hrwsi.py --dry-run   # list, don't fetch
    uv run python scripts/1b-download_data-hrwsi.py

Always do a --dry-run first: it reports how many products and how many GB the
query matches, and prints the layer filenames actually present so that
`hrwsi.layer_patterns` can be narrowed before committing to a bulk download.

Output mirrors the S3 key layout so re-runs are resumable and provenance is kept:
    derived/eo/hrwsi/<PRODUCT>/<tile>/<YYYY>/<MM>/<DD>/<product_name>/<layers>
plus a manifest CSV listing every file retrieved.

Legal notice: Copernicus data are free, full and open (Regulation (EU) No
1159/2013). Publications must state the source and that the data were produced
with funding by the European Union, and must flag any modification.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import os
import sqlite3
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import date, timedelta
from pathlib import Path

import yaml

repo_root = Path(__file__).resolve().parent.parent
SPEC_PATH = repo_root / "inputs" / "data_specifications.yaml"

# Public HR-WSI dissemination bucket (see module docstring).
DEFAULT_ENDPOINT = "https://s3.WAW3-2.cloudferro.com"
DEFAULT_BUCKET = "HRWSI"
DEFAULT_ACCESS_KEY = "c4ae60af7b144053803c618a8860f7c9"
DEFAULT_SECRET_KEY = "dcb3ba1f6eab45aaaec5802feef5e2e4"

# Daily-cadence product types. The yearly ones (SP_S2, ICD, WCD, ...) are keyed by
# hydrological year instead of Y/M/D and are out of scope here.
DAILY_PRODUCTS = {"FSC", "SWS", "GFSC", "WDS", "CC", "WIC_S1", "WIC_S2", "WIC_S1S2"}

# HR-WSI starts here; earlier dates return nothing.
ARCHIVE_START = date(2016, 9, 1)

_DEFAULTS = {
    "product_types": ["SWS", "GFSC"],
    "tiles": [],
    "mgrs_gpkg": "inputs/MGRS_tiles.gpkg",
    "first_date": "2016-09-01",
    "last_date": None,  # None -> today
    "months": [10, 11, 12, 1, 2, 3, 4, 5, 6],
    "layer_patterns": ["*"],
    "output_dir": "derived/eo/hrwsi",
    "max_workers": 4,
    "endpoint_url": DEFAULT_ENDPOINT,
    "bucket": DEFAULT_BUCKET,
}


def load_config() -> dict:
    """Read inputs/data_specifications.yaml -> hrwsi block; fall back to defaults."""
    cfg = dict(_DEFAULTS)
    if not SPEC_PATH.is_file():
        print(f"Note: {SPEC_PATH} missing; using built-in defaults.", flush=True)
        return cfg
    with SPEC_PATH.open(encoding="utf-8") as f:
        spec = yaml.safe_load(f) or {}
    for key, value in (spec.get("hrwsi") or {}).items():
        if value is not None:
            cfg[key] = value
    # The bbox is shared with the ERA5-Land download so the two stay on one domain.
    cfg["bbox"] = (spec.get("download") or {}).get("bbox")
    return cfg


def parse_date(text: str | date) -> date:
    if isinstance(text, date):
        return text
    return date.fromisoformat(str(text))


def tiles_from_bbox(bbox: list[float], gpkg_path: Path) -> list[str]:
    """Look up the Sentinel-2 / HR-WSI MGRS tiles intersecting a lon/lat bbox.

    Reads the MGRS_tiles.gpkg shipped with the official EEA client using the
    GeoPackage R-tree index, so only the standard library is needed (no
    geopandas). The R-tree holds bounding boxes, so the result can be very
    slightly over-inclusive -- harmless, since a tile with no data simply
    returns nothing from S3.
    """
    if not gpkg_path.is_file():
        raise SystemExit(
            f"Cannot derive tiles: {gpkg_path} not found.\n"
            "Either set hrwsi.tiles explicitly in inputs/data_specifications.yaml, "
            "or download MGRS_tiles.gpkg from\n"
            "  https://github.com/eea/clms-hrwsi-api-client-python\n"
            f"and place it at {gpkg_path}."
        )
    xmin, ymin, xmax, ymax = bbox
    con = sqlite3.connect(f"file:{gpkg_path}?mode=ro", uri=True)
    try:
        rows = con.execute(
            """
            SELECT t.Name
              FROM sentinel_2_index_shapefile t
              JOIN rtree_sentinel_2_index_shapefile_geom r ON t.fid = r.id
             WHERE r.maxx >= ? AND r.minx <= ? AND r.maxy >= ? AND r.miny <= ?
            """,
            (xmin, xmax, ymin, ymax),
        ).fetchall()
    finally:
        con.close()
    return sorted({str(row[0]) for row in rows})


def month_starts(first: date, last: date) -> list[tuple[int, int]]:
    """Every (year, month) touched by the date range, in order."""
    out = []
    year, month = first.year, first.month
    while (year, month) <= (last.year, last.month):
        out.append((year, month))
        month += 1
        if month == 13:
            year, month = year + 1, 1
    return out


def make_bucket(cfg: dict):
    try:
        import boto3
    except ImportError:
        raise SystemExit(
            "boto3 is required. Add it with `uv add boto3` (it is already declared "
            "in pyproject.toml, so `uv sync` should be enough)."
        )
    session = boto3.resource(
        "s3",
        endpoint_url=cfg["endpoint_url"],
        aws_access_key_id=os.environ.get("HRWSI_S3_ACCESS_KEY", DEFAULT_ACCESS_KEY),
        aws_secret_access_key=os.environ.get("HRWSI_S3_SECRET_KEY", DEFAULT_SECRET_KEY),
    )
    return session.Bucket(cfg["bucket"])


def list_month(bucket, product_type: str, tile: str, year: int, month: int) -> list:
    """All object summaries under one product/tile/year/month prefix."""
    prefix = f"{product_type}/{tile}/{year}/{month:02d}/"
    return list(bucket.objects.filter(Prefix=prefix))


def day_of_key(key: str) -> date | None:
    """Parse the YYYY/MM/DD embedded in an HR-WSI S3 key."""
    parts = key.split("/")
    if len(parts) < 5:
        return None
    try:
        return date(int(parts[2]), int(parts[3]), int(parts[4]))
    except (ValueError, IndexError):
        return None


def matches_layer(key: str, patterns: list[str]) -> bool:
    name = key.rsplit("/", 1)[-1]
    return any(fnmatch.fnmatch(name, pattern) for pattern in patterns)


def build_query(bucket, cfg: dict, tiles: list[str], first: date, last: date) -> list[dict]:
    """List every matching object across product types, tiles and months."""
    months = set(cfg["months"]) if cfg["months"] else set(range(1, 13))
    patterns = cfg["layer_patterns"] or ["*"]
    wanted_months = [(y, m) for (y, m) in month_starts(first, last) if m in months]

    jobs = [
        (pt, tile, y, m)
        for pt in cfg["product_types"]
        for tile in tiles
        for (y, m) in wanted_months
    ]
    print(
        f"Scanning {len(jobs)} product/tile/month prefixes "
        f"({len(cfg['product_types'])} product types x {len(tiles)} tiles x "
        f"{len(wanted_months)} months)",
        flush=True,
    )

    found: list[dict] = []
    with ThreadPoolExecutor(max_workers=int(cfg["max_workers"])) as pool:
        futures = {
            pool.submit(list_month, bucket, pt, tile, y, m): (pt, tile, y, m)
            for (pt, tile, y, m) in jobs
        }
        for done, future in enumerate(as_completed(futures), start=1):
            pt, tile, y, m = futures[future]
            try:
                objects = future.result()
            except Exception as err:  # noqa: BLE001 - report and carry on
                print(f"  ! {pt}/{tile}/{y}/{m:02d}: {err}", flush=True)
                continue
            for obj in objects:
                day = day_of_key(obj.key)
                if day is None or not (first <= day <= last):
                    continue
                if not matches_layer(obj.key, patterns):
                    continue
                found.append({"key": obj.key, "size": int(obj.size)})
            if done % 50 == 0 or done == len(jobs):
                print(f"  scanned {done}/{len(jobs)} prefixes, {len(found)} files so far", flush=True)
    return found


def summarise(found: list[dict]) -> None:
    products = sorted({item["key"].rsplit("/", 1)[0] for item in found})
    total_gb = sum(item["size"] for item in found) / 1e9
    print(f"\nMatched {len(found)} files in {len(products)} products, {total_gb:.2f} GB")
    if not found:
        return
    layers: dict[str, int] = {}
    for item in found:
        layers[item["key"].rsplit("/", 1)[-1].split("_")[-1]] = (
            layers.get(item["key"].rsplit("/", 1)[-1].split("_")[-1], 0) + 1
        )
    print("Layer file suffixes present (use these to narrow hrwsi.layer_patterns):")
    for suffix, count in sorted(layers.items(), key=lambda kv: -kv[1]):
        print(f"  {suffix:<28} {count}")
    print("\nFirst few products:")
    for product in products[:5]:
        print(f"  {product}")


def download_one(bucket, key: str, size: int, out_root: Path) -> tuple[str, str]:
    """Fetch one object, skipping it if a same-size file is already on disk."""
    target = out_root / key
    if target.is_file() and target.stat().st_size == size:
        return key, "skipped"
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    bucket.download_file(key, str(partial))
    partial.replace(target)
    return key, "downloaded"


def main() -> None:
    cfg = load_config()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true", help="list matching products and stop")
    parser.add_argument("--product-types", nargs="+", help="override hrwsi.product_types")
    parser.add_argument("--tiles", nargs="+", help="override hrwsi.tiles (e.g. 32TPS 32TPT)")
    parser.add_argument("--date-start", help="override hrwsi.first_date (YYYY-MM-DD)")
    parser.add_argument("--date-end", help="override hrwsi.last_date (YYYY-MM-DD)")
    parser.add_argument("--months", nargs="+", type=int, help="override hrwsi.months")
    parser.add_argument("--output-dir", help="override hrwsi.output_dir")
    args = parser.parse_args()

    if args.product_types:
        cfg["product_types"] = args.product_types
    if args.tiles:
        cfg["tiles"] = args.tiles
    if args.date_start:
        cfg["first_date"] = args.date_start
    if args.date_end:
        cfg["last_date"] = args.date_end
    if args.months:
        cfg["months"] = args.months
    if args.output_dir:
        cfg["output_dir"] = args.output_dir

    unknown = [pt for pt in cfg["product_types"] if pt not in DAILY_PRODUCTS]
    if unknown:
        raise SystemExit(
            f"Unsupported product type(s): {unknown}. This script handles the daily "
            f"products {sorted(DAILY_PRODUCTS)}; the yearly ones (SP_S2, ICD, WCD, ...) "
            "are keyed by hydrological year and need different prefixes."
        )

    first = parse_date(cfg["first_date"])
    last = parse_date(cfg["last_date"]) if cfg["last_date"] else date.today()
    if first < ARCHIVE_START:
        print(f"Note: HR-WSI starts {ARCHIVE_START}; clamping start date.", flush=True)
        first = ARCHIVE_START
    if first > last:
        raise SystemExit(f"first_date {first} is after last_date {last}")

    tiles = list(cfg["tiles"] or [])
    if not tiles:
        if not cfg.get("bbox"):
            raise SystemExit("No hrwsi.tiles set and no download.bbox to derive them from.")
        tiles = tiles_from_bbox(cfg["bbox"], repo_root / cfg["mgrs_gpkg"])
        print(f"Derived {len(tiles)} tiles from bbox {cfg['bbox']}: {' '.join(tiles)}", flush=True)
    tiles = [tile[1:] if len(tile) == 6 and tile.upper().startswith("T") else tile for tile in tiles]

    print(f"Products : {' '.join(cfg['product_types'])}")
    print(f"Tiles    : {' '.join(tiles)}")
    print(f"Dates    : {first} .. {last}  (months {sorted(set(cfg['months']))})")
    print(f"Layers   : {cfg['layer_patterns']}")

    bucket = make_bucket(cfg)
    found = build_query(bucket, cfg, tiles, first, last)
    summarise(found)

    if args.dry_run:
        print("\n--dry-run: nothing downloaded.")
        return
    if not found:
        return

    out_root = repo_root / cfg["output_dir"]
    out_root.mkdir(parents=True, exist_ok=True)
    print(f"\nDownloading into {out_root}", flush=True)

    results: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=int(cfg["max_workers"])) as pool:
        futures = {
            pool.submit(download_one, bucket, item["key"], item["size"], out_root): item
            for item in found
        }
        for done, future in enumerate(as_completed(futures), start=1):
            item = futures[future]
            try:
                results.append(future.result())
            except Exception as err:  # noqa: BLE001 - one bad object must not sink the run
                print(f"  ! {item['key']}: {err}", flush=True)
                results.append((item["key"], f"failed: {err}"))
            if done % 100 == 0 or done == len(found):
                print(f"  {done}/{len(found)} files", flush=True)

    manifest = out_root / "manifest.csv"
    write_header = not manifest.is_file()
    with manifest.open("a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["key", "product_type", "tile", "date", "status"])
        for key, status in sorted(results):
            parts = key.split("/")
            day = day_of_key(key)
            writer.writerow([key, parts[0], parts[1], day.isoformat() if day else "", status])

    counts: dict[str, int] = {}
    for _, status in results:
        label = status.split(":")[0]
        counts[label] = counts.get(label, 0) + 1
    print(f"\nDone: {counts}. Manifest: {manifest}")


if __name__ == "__main__":
    main()
