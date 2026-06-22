"""
save as download_viirs_deep_blue.py
Requires: pip install requests python-cmr tqdm
Optional for convenience: pip install earthaccess
"""

import os
import netrc
import requests
from datetime import datetime
from urllib.parse import urlsplit
from tqdm import tqdm

CMR_GRANULES_UMM = "https://cmr.earthdata.nasa.gov/search/granules.umm_json"

# Earthdata recommended session with header redirection (from EDL docs)
class SessionWithHeaderRedirection(requests.Session):
    AUTH_HOST = "urs.earthdata.nasa.gov"

    def __init__(self, username, password):
        super().__init__()
        self.auth = (username, password)

    def rebuild_auth(self, prepared_request, response):
        headers = prepared_request.headers
        url = prepared_request.url
        # Remove Authorization header when crossing hosts unless the host is UAS auth host
        if "Authorization" in headers:
            original_parsed = requests.utils.urlparse(response.request.url)
            redirect_parsed = requests.utils.urlparse(url)
            if (original_parsed.hostname != redirect_parsed.hostname) and \
               (redirect_parsed.hostname != self.AUTH_HOST) and \
               (original_parsed.hostname != self.AUTH_HOST):
                del headers["Authorization"]
        return

def _load_netrc_credentials(host="urs.earthdata.nasa.gov"):
    try:
        auth = netrc.netrc()
        creds = auth.authenticators(host)
        if creds:
            username, _, password = creds
            return username, password
    except Exception:
        return None, None
    return None, None

def _page_cmr(short_name, west, south, east, north, date_iso, page_size=200):
    # returns generator of page JSON 'items'
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
        for it in items:
            yield it
        # if returned less than page_size, we are done
        if len(items) < page_size:
            break
        page += 1

def _choose_download_url(related_urls, prefer_nc=True):
    """
    Pick a download URL from RelatedUrls list.
    - prefer_nc=True: prefer links ending in .nc
    Preference order:
      1) HTTPS .nc links
      2) HTTPS any other links
      3) S3 .nc links
      4) Other S3 or OPENDAP links
    """
    https_nc = []
    https_any = []
    s3_nc = []
    s3_any = []
    opendap = []

    for r in related_urls:
        url = r.get("URL") or r.get("url")
        if not url:
            continue
        typ = (r.get("Type") or "").upper()
        subtype = (r.get("Subtype") or "").upper()

        if url.startswith("https://"):
            if url.endswith(".nc"):
                https_nc.append(url)
            else:
                https_any.append(url)
        elif url.startswith("s3://"):
            if url.endswith(".nc"):
                s3_nc.append(url)
            else:
                s3_any.append(url)
        elif "OPENDAP" in subtype or "OPENDAP" in typ or "USE SERVICE" in typ:
            opendap.append(url)

    if prefer_nc:
        if https_nc:
            return https_nc[0]
        if s3_nc:
            return s3_nc[0]
    # fall back if no .nc links
    if https_any:
        return https_any[0]
    if s3_any:
        return s3_any[0]
    if opendap:
        return opendap[0]
    return None


def download_viirs_deep_blue(bbox, date, out_dir="downloads", username=None, password=None,
                             short_name="AERDB_L2_VIIRS_SNPP", verbose=True):
    """
    bbox: (west, south, east, north) floats in degrees
    date: string 'YYYY-MM-DD' (the function searches that UTC calendar day)
    out_dir: directory to save files
    username/password: optional; if None the function tries ~/.netrc for EDL creds
    short_name: collection short name for VIIRS SNPP Deep Blue (default is AERDB_L2_VIIRS_SNPP)
    Returns: list of downloaded file paths
    """
    os.makedirs(out_dir, exist_ok=True)

    # load credentials if not provided
    if (username is None) or (password is None):
        u, p = _load_netrc_credentials()
        if u and p:
            username, password = username or u, password or p
    if not (username and password):
        raise RuntimeError("No Earthdata credentials found: either pass username/password or set ~/.netrc for urs.earthdata.nasa.gov")

    west, south, east, north = bbox
    # quick date validation
    try:
        datetime.fromisoformat(date)
    except Exception:
        raise ValueError("date must be YYYY-MM-DD")

    downloaded = []
    session = SessionWithHeaderRedirection(username, password)

    # iterate CMR granule pages
    if verbose:
        print(f"Querying CMR for {short_name} intersecting bbox {bbox} on {date}...")

    for gran_item in _page_cmr(short_name, west, south, east, north, date_iso=date):
        # UMM-G JSON layout: item['umm'] contains granule-level info including RelatedUrls
        umm = gran_item.get("umm") or {}
        related = umm.get("RelatedUrls") or []
        if not related:
            # try legacy structure
            # some providers used different keys; try the 'meta' or 'granule' fallback
            # continue if no links
            continue

        url = _choose_download_url(related)
        if not url:
            # nothing to download
            continue

        # derive filename
        parsed = urlsplit(url)
        filename = os.path.basename(parsed.path) or parsed.path.replace("/", "_")
        local_path = os.path.join(out_dir, filename)

        # Skip if already downloaded
        if os.path.exists(local_path):
            if verbose:
                print(f"Skipping (exists): {filename}")
            downloaded.append(local_path)
            continue

        # If s3:// - warn (not directly handled here)
        if url.startswith("s3://"):
            # either use boto3 (if you have AWS credentials) or `aws s3 cp` externally
            print(f"Found S3 URL for granule: {url}. This function currently prefers https downloads. Use boto3 or awscli to fetch S3 objects.")
            continue

        # download with streaming & auth
        if verbose:
            print(f"Downloading: {url} -> {local_path}")
        try:
            with session.get(url, stream=True, allow_redirects=True, timeout=120) as r:
                r.raise_for_status()
                total = int(r.headers.get("Content-Length", 0))
                chunk = 1024 * 1024
                with open(local_path, "wb") as f, tqdm(total=total, unit="B", unit_scale=True, desc=filename, disable=not verbose) as pbar:
                    for chunk_bytes in r.iter_content(chunk_size=chunk):
                        if chunk_bytes:
                            f.write(chunk_bytes)
                            pbar.update(len(chunk_bytes))
            downloaded.append(local_path)
        except Exception as e:
            print(f"Failed to download {url}: {e}")

    if verbose:
        print("Done. Downloaded files:", downloaded)
    return downloaded

# Example usage:
if __name__ == "__main__":
    import argparse, json
    parser = argparse.ArgumentParser(description="Download VIIRS Deep Blue Aerosol data")
    parser.add_argument("--bbox", nargs=4, type=float, required=True,
                        help="Bounding box: west south east north")
    parser.add_argument("--date", type=str, required=True,
                        help="Date (YYYY-MM-DD)")
    parser.add_argument("--outdir", type=str, default="downloads",
                        help="Directory to save files")
    args = parser.parse_args()

    files = download_viirs_deep_blue(
        bbox=(args.bbox[0], args.bbox[1], args.bbox[2], args.bbox[3]),
        date=args.date,
        out_dir=args.outdir
    )
    # print list of downloaded file paths as JSON
    print(json.dumps(files))
