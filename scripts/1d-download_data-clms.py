"""
Download CLMS bio-geophysical products (soil moisture, snow cover) from CDSE OData.

Third of the Earth-observation downloaders, alongside 1b (HR-WSI snow rasters) and
1c (Sentinel-1 wet-snow fraction). These are the *other* preconditioning fields --
how wet the ground already was, and daily snow cover with better temporal
continuity than the Sentinel-2 products:

  swi   Soil Water Index, 1 km, Europe, daily. Root-zone soil moisture from
        SCATSAR (Sentinel-1 + ASCAT). The antecedent-wetness variable: it carries
        past observations forward in time, so it degrades more gracefully than raw
        surface soil moisture.
  ssm   Surface Soil Moisture, 1 km, Europe, daily, Sentinel-1.
  sce   Snow Cover Extent, 1 km, global, daily, Sentinel-3 SLSTR (ENVEO).

Honest limitation for the soil moisture pair: retrieval is not possible over
snow-covered or frozen ground, so these fields are masked during exactly the
conditions a rain-on-snow event happens in. Their value here is the *antecedent*
state -- autumn and early-winter wetness before the pack establishes, and the
shoulder-season events -- not the state during the event.

Unlike 1b, this one needs a free Copernicus Data Space Ecosystem account
(https://dataspace.copernicus.eu). Supply credentials via environment variables:

    export CDSE_USERNAME=you@example.org
    export CDSE_PASSWORD='...'

Usage:
    uv run python scripts/1d-download_data-clms.py --dry-run
    uv run python scripts/1d-download_data-clms.py

Run --dry-run first. These are pan-European / global daily NetCDF files at 1 km,
so a full archive is large; narrow `clms.months` and the date range before
committing.

Output:
    derived/eo/clms/<dataset>/<YYYY>/<filename>.nc

Legal notice: Copernicus data are free, full and open (Regulation (EU) No
1159/2013). Publications must state the source and that the data were produced
with funding by the European Union, and must flag any modification.
"""

from __future__ import annotations

import argparse
import json
import os
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

import yaml

repo_root = Path(__file__).resolve().parent.parent
SPEC_PATH = repo_root / "inputs" / "data_specifications.yaml"

TOKEN_URL = (
    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/"
    "protocol/openid-connect/token"
)
CATALOGUE_URL = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
DOWNLOAD_URL = "https://download.dataspace.copernicus.eu/odata/v1/Products({id})/$value"

# datasetAlias values as they appear in the CDSE CLMS catalogue.
DATASETS = {
    "swi": "swi_europe_1km_daily",
    "ssm": "ssm_europe_1km_daily",
    "sce": "sce_global_1km_daily",
}

_DEFAULTS = {
    "datasets": ["swi", "sce"],
    "first_date": "2016-09-01",
    "last_date": None,  # None -> today
    "months": [10, 11, 12, 1, 2, 3, 4, 5, 6],
    "output_dir": "derived/eo/clms",
    "page_size": 500,
}


def load_config() -> dict:
    cfg = dict(_DEFAULTS)
    if not SPEC_PATH.is_file():
        print(f"Note: {SPEC_PATH} missing; using built-in defaults.", flush=True)
        return cfg
    with SPEC_PATH.open(encoding="utf-8") as f:
        spec = yaml.safe_load(f) or {}
    for key, value in (spec.get("clms") or {}).items():
        if value is not None:
            cfg[key] = value
    return cfg


def get_token() -> str:
    username = os.environ.get("CDSE_USERNAME")
    password = os.environ.get("CDSE_PASSWORD")
    if not username or not password:
        raise SystemExit(
            "Set CDSE_USERNAME and CDSE_PASSWORD (free account at "
            "https://dataspace.copernicus.eu). The catalogue is open, but "
            "downloading a product needs a token."
        )
    payload = urllib.parse.urlencode(
        {
            "grant_type": "password",
            "username": username,
            "password": password,
            "client_id": "cdse-public",
        }
    ).encode()
    request = urllib.request.Request(TOKEN_URL, data=payload)
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.load(response)["access_token"]
    except urllib.error.HTTPError as err:
        raise SystemExit(f"CDSE token request failed ({err.code}): {err.read().decode()[:400]}")


def query_products(alias: str, first: date, last: date, page_size: int) -> list[dict]:
    """Every CLMS product of one dataset within the date range (paged)."""
    odata_filter = (
        "Collection/Name eq 'CLMS' and "
        "Attributes/OData.CSC.StringAttribute/any("
        f"att:att/Name eq 'datasetAlias' and att/OData.CSC.StringAttribute/Value eq '{alias}') and "
        f"ContentDate/Start ge {first.isoformat()}T00:00:00.000Z and "
        f"ContentDate/Start lt {last.isoformat()}T23:59:59.999Z"
    )
    url = (
        CATALOGUE_URL
        + "?"
        + urllib.parse.urlencode(
            {"$filter": odata_filter, "$top": page_size, "$orderby": "ContentDate/Start asc"}
        )
    )

    products: list[dict] = []
    while url:
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                payload = json.load(response)
        except urllib.error.HTTPError as err:
            raise SystemExit(
                f"CDSE catalogue query failed ({err.code}) for {alias}: "
                f"{err.read().decode()[:400]}"
            )
        products.extend(payload.get("value", []))
        url = payload.get("@odata.nextLink")
        print(f"  {alias}: {len(products)} products listed", flush=True)
    return products


def product_date(product: dict) -> date | None:
    stamp = (product.get("ContentDate") or {}).get("Start")
    if not stamp:
        return None
    try:
        return date.fromisoformat(stamp[:10])
    except ValueError:
        return None


def download_product(product: dict, token: str, out_dir: Path) -> str:
    name = product["Name"]
    if not name.endswith(".nc"):
        name = f"{name}.nc"
    target = out_dir / name
    size = int(product.get("ContentLength") or 0)
    if target.is_file() and (size == 0 or target.stat().st_size == size):
        return "skipped"

    url = DOWNLOAD_URL.format(id=product["Id"])
    request = urllib.request.Request(url, headers={"Authorization": f"Bearer {token}"})
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    with urllib.request.urlopen(request, timeout=600) as response, partial.open("wb") as handle:
        while True:
            block = response.read(1 << 20)
            if not block:
                break
            handle.write(block)
    partial.replace(target)
    return "downloaded"


def main() -> None:
    cfg = load_config()
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dry-run", action="store_true", help="list matching products and stop")
    parser.add_argument("--datasets", nargs="+", choices=sorted(DATASETS), help="override clms.datasets")
    parser.add_argument("--date-start", help="override clms.first_date (YYYY-MM-DD)")
    parser.add_argument("--date-end", help="override clms.last_date (YYYY-MM-DD)")
    parser.add_argument("--months", nargs="+", type=int, help="override clms.months")
    args = parser.parse_args()

    if args.datasets:
        cfg["datasets"] = args.datasets
    if args.date_start:
        cfg["first_date"] = args.date_start
    if args.date_end:
        cfg["last_date"] = args.date_end
    if args.months:
        cfg["months"] = args.months

    unknown = [name for name in cfg["datasets"] if name not in DATASETS]
    if unknown:
        raise SystemExit(f"Unknown dataset(s) {unknown}; choose from {sorted(DATASETS)}.")

    first = date.fromisoformat(str(cfg["first_date"]))
    last = date.fromisoformat(str(cfg["last_date"])) if cfg["last_date"] else date.today()
    months = set(cfg["months"]) if cfg["months"] else set(range(1, 13))
    out_root = repo_root / cfg["output_dir"]

    selected: dict[str, list[dict]] = {}
    for name in cfg["datasets"]:
        print(f"Querying {name} ({DATASETS[name]}) {first} .. {last}", flush=True)
        products = query_products(DATASETS[name], first, last, int(cfg["page_size"]))
        keep = [p for p in products if (d := product_date(p)) and d.month in months]
        total_gb = sum(int(p.get("ContentLength") or 0) for p in keep) / 1e9
        print(f"  {name}: {len(keep)} products in the requested months, {total_gb:.2f} GB")
        selected[name] = keep

    if args.dry_run:
        print("\n--dry-run: nothing downloaded.")
        for name, products in selected.items():
            for product in products[:3]:
                print(f"  {name}: {product['Name']}")
        return

    token = get_token()
    for name, products in selected.items():
        counts: dict[str, int] = {}
        for done, product in enumerate(products, start=1):
            day = product_date(product)
            out_dir = out_root / name / (str(day.year) if day else "unknown")
            try:
                status = download_product(product, token, out_dir)
            except urllib.error.HTTPError as err:
                if err.code in (401, 403):
                    # Tokens are short lived; refresh once and retry this product.
                    token = get_token()
                    status = download_product(product, token, out_dir)
                else:
                    print(f"  ! {product['Name']}: HTTP {err.code}", flush=True)
                    status = "failed"
            except Exception as err:  # noqa: BLE001 - one bad file must not sink the run
                print(f"  ! {product['Name']}: {err}", flush=True)
                status = "failed"
            counts[status] = counts.get(status, 0) + 1
            if done % 25 == 0 or done == len(products):
                print(f"  {name}: {done}/{len(products)} {counts}", flush=True)
        print(f"{name}: {counts}")


if __name__ == "__main__":
    main()
