#!/usr/bin/env python3

"""
VIIRS Deep Blue downloader
Returns JSON list of downloaded files to stdout.
Requires: requests, tqdm, netrc
"""

import os
import sys
import json
import netrc
import requests
from datetime import datetime
from urllib.parse import urlsplit

try:
    from tqdm import tqdm
except Exception:
    class _Dummy:
        def __init__(self, *a, **k): pass
        def update(self, n): pass
        def __enter__(self): return self
        def __exit__(self, *a): pass
    def tqdm(*a, **k): return _Dummy()

CMR_GRANULES_UMM = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"


# ----------------------------------------------------------------------
# AUTH SESSION — USES ~/.netrc (NO HARDCODED PASSWORDS)
# ----------------------------------------------------------------------
class SessionWithHeaderRedirection(requests.Session):
    AUTH_HOST = "urs.earthdata.nasa.gov"

    def __init__(self):
        super().__init__()
        self.auth = self._load_credentials()

    def _load_credentials(self):
        try:
            creds = netrc.netrc().authenticators(self.AUTH_HOST)
            if creds:
                username, _, password = creds
                return (username, password)
        except Exception:
            pass
        raise RuntimeError(
            "No Earthdata credentials found. "
            "Create ~/.netrc with:\n"
            "machine urs.earthdata.nasa.gov login YOUR_USER password YOUR_PASS"
        )

    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        orig = requests.utils.urlparse(response.request.url)
        redir = requests.utils.urlparse(prepared_request.url)

        if ("Authorization" in headers and
            orig.hostname != redir.hostname and
            redir.hostname != self.AUTH_HOST and
            orig.hostname != self.AUTH_HOST):
            del headers["Authorization"]
        return


# ----------------------------------------------------------------------
# Helper: CMR paging
# ----------------------------------------------------------------------
def cmr_search(short_name, bbox, date_iso, page_size=200):
    west, south, east, north = bbox
    page = 1

    while True:
        params = {
            "short_name": short_name,
            "bounding_box": f"{west},{south},{east},{north}",
            "temporal": f"{date_iso}T00:00:00Z,{date_iso}T23:59:59Z",
            "page_size": page_size,
            "page_num": page,
        }

        r = requests.get(CMR_GRANULES_UMM, params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        items = data.get("items", [])

        if not items:
            break

        for item in items:
            yield item

        if len(items) < page_size:
            break

        page += 1


# ----------------------------------------------------------------------
# Helper: choose best URL
# ----------------------------------------------------------------------
def choose_download_url(related_urls):
    https_nc = []
    https_any = []

    for r in related_urls:
        url = r.get("URL") or r.get("url")
        if not url:
            continue
        if url.startswith("https://"):
            if url.endswith(".nc") or ".nc4" in url:
                https_nc.append(url)
            else:
                https_any.append(url)

    if https_nc:
        return https_nc[0]
    if https_any:
        return https_any[0]
    return None


# ----------------------------------------------------------------------
# Downloader core
# ----------------------------------------------------------------------
def download_viirs(bbox, date, outdir, platforms):
    os.makedirs(outdir, exist_ok=True)

    platforms = [p.strip().upper() for p in platforms.split(",")]
    platform_map = {
        "SNPP": "AERDB_L2_VIIRS_SNPP",
        "NOAA20": "AERDB_L2_VIIRS_NOAA20",
        "NOAA21": "AERDB_L2_VIIRS_NOAA21",
    }
    short_names = [platform_map[p] for p in platforms if p in platform_map]

    if not short_names:
        raise ValueError("No valid platforms specified (SNPP, NOAA20, NOAA21).")

    try:
        datetime.fromisoformat(date)
    except Exception:
        raise ValueError("Date must be YYYY-MM-DD")

    session = SessionWithHeaderRedirection()
    downloaded = []

    for short in short_names:
        print(f"Searching: {short}", file=sys.stderr)
        found_any = False

        for gran in cmr_search(short, bbox, date):
            found_any = True
            umm = gran.get("umm", {})
            urls = umm.get("RelatedUrls", [])
            url = choose_download_url(urls)

            if not url:
                continue

            fname = os.path.basename(urlsplit(url).path)
            local_path = os.path.join(outdir, fname)

            if os.path.exists(local_path):
                downloaded.append(local_path)
                continue

            print(f"Downloading {url} -> {local_path}", file=sys.stderr)
            try:
                with session.get(url, stream=True, timeout=120) as r:
                    r.raise_for_status()
                    total = int(r.headers.get("Content-Length", 0))
                    chunk = 1024 * 1024
                    with open(local_path, "wb") as f, tqdm(total=total, unit="B", unit_scale=True) as pbar:
                        for chunk_bytes in r.iter_content(chunk_size=chunk):
                            if chunk_bytes:
                                f.write(chunk_bytes)
                                pbar.update(len(chunk_bytes))
                downloaded.append(local_path)
            except Exception as e:
                print(f"Failed to download {url}: {e}", file=sys.stderr)

        if not found_any:
            print(f"No granules found for {short}", file=sys.stderr)

    # Return JSON ALWAYS
    print(json.dumps(downloaded))
    return downloaded


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--bbox", nargs=4, type=float, required=True,
                        help="Bounding box: W S E N")
    parser.add_argument("--date", type=str, required=True)
    parser.add_argument("--outdir", type=str, default="downloads")
    parser.add_argument("--platforms", type=str, default="NOAA20,NOAA21")

    args = parser.parse_args()

    download_viirs(
        bbox=(args.bbox[0], args.bbox[1], args.bbox[2], args.bbox[3]),
        date=args.date,
        outdir=args.outdir,
        platforms=args.platforms
    )
