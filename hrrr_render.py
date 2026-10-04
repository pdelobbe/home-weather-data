"""Render HRRR overlays (reflectivity+ptype, STP, UH 2-5 km run-max) for central Oklahoma.

Downloads only the needed GRIB2 messages (HTTP byte ranges via the .idx files) from
NOAA's HRRR bucket on AWS, regrids them bilinearly onto a Web Mercator image grid,
and writes transparent PNGs + meta.json for Leaflet imageOverlay.

Usage: python hrrr_render.py OUT_DIR [--hours 18] [--run YYYYMMDDHH]

OUT_DIR/{refl,stp,uh} are cleared each run (only the current run is kept); meta.json is
written last, atomically. A forecast hour that fails to download is retried twice; the
run then stops at the last good hour (meta.hours), failing if fewer than 6 hours succeed.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import eccodes
import numpy as np
from PIL import Image
from pyproj import Proj
from scipy.ndimage import gaussian_filter, map_coordinates

BUCKET = "https://noaa-hrrr-bdp-pds.s3.amazonaws.com"
# Bounds of the rendered image (lon/lat). Centered near home (35.615, -97.693).
WEST, EAST, SOUTH, NORTH = -100.2, -95.2, 33.9, 37.3
WIDTH = 1200  # px; ~370 m per pixel
MIN_HOURS = 6  # fail the run if fewer complete forecast hours than this (or --hours, if smaller)

FIELDS = {
    "refd": "REFD:1000 m above ground",
    "crain": "CRAIN:surface",
    "csnow": "CSNOW:surface",
    "cfrzr": "CFRZR:surface",
    "cicep": "CICEP:surface",
    "uh25": "MXUPHL:5000-2000 m above ground",
    "cape": "CAPE:surface",
    "cin": "CIN:surface",
    "srh1": "HLCY:1000-0 m above ground",
    "ushr6": "VUCSH:0-6000 m above ground",
    "vshr6": "VVCSH:0-6000 m above ground",
    "t2": "TMP:2 m above ground",
    "td2": "DPT:2 m above ground",
}

# (threshold, hex) — value >= threshold gets that color
RAIN = [(10, "#a8f0a0"), (20, "#4cc94a"), (30, "#1e9a2e"), (35, "#f7e83a"), (40, "#f5b62a"),
        (45, "#f07c22"), (50, "#e3301f"), (55, "#b5121b"), (60, "#d63fd0"), (65, "#8b4fd6")]
SNOW = [(10, "#b9dcff"), (20, "#7fb6f5"), (30, "#3f86e3"), (40, "#1f4fb8")]
FRZR = [(10, "#f9c0d8"), (20, "#f08fbd"), (30, "#de4f9a")]
SLEET = [(10, "#d8bff2"), (20, "#b48ae6"), (30, "#8a52cf")]
# Continuous run-max UH scale (m2/s2): weak rotation recedes in grays so supercell
# swaths stand out — blue 50-100, yellow→orange 100-150, red 150-200, purple 200-300, pink 300+.
UH_STOPS = [(25, "#bdbdbd"), (50, "#7d7d7d"), (50.01, "#34587a"), (100, "#9cc7cf"),
            (100.01, "#f3df8f"), (150, "#e8700f"), (150.01, "#e0451b"), (200, "#8d0f45"),
            (200.01, "#6f1d9c"), (300, "#d58ae8"), (300.01, "#f2bccb"), (400, "#b5385a")]
UH_ALPHA = [(25, 60), (50, 150), (100, 240)]
STP = [(0.5, "#bfe3ff"), (1, "#6fb3f0"), (2, "#f2d43a"), (3, "#f39a2c"), (4, "#e2412a"),
       (6, "#b51b52"), (8, "#c44fc4")]


def hexrgb(h: str) -> tuple[int, int, int]:
    return int(h[1:3], 16), int(h[3:5], 16), int(h[5:7], 16)


def with_retries(fn, *args, retries: int = 2, backoff: float = 2.0):
    """Call fn(*args); on network/HTTP errors retry `retries` times with exponential backoff."""
    for attempt in range(retries + 1):
        try:
            return fn(*args)
        except (urllib.error.URLError, OSError):  # HTTPError is a URLError; timeouts are OSError
            if attempt == retries:
                raise
            time.sleep(backoff * 2 ** attempt)


def fetch(url: str, rng: str | None = None) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "HomeWeather/1.0"})
    if rng:
        req.add_header("Range", f"bytes={rng}")
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def url_for(run: dt.datetime, fh: int, ext: str = "") -> str:
    return f"{BUCKET}/hrrr.{run:%Y%m%d}/conus/hrrr.t{run:%H}z.wrfsfcf{fh:02d}.grib2{ext}"


def latest_run(hours: int) -> dt.datetime:
    now = dt.datetime.now(dt.timezone.utc).replace(minute=0, second=0, microsecond=0, tzinfo=None)
    for back in range(1, 6):
        run = now - dt.timedelta(hours=back)
        try:
            fetch(url_for(run, hours, ".idx"))
            return run
        except Exception:
            continue
    raise SystemExit("no complete HRRR run found in the last 5 hours")


def byte_ranges(idx: str) -> dict[str, str]:
    lines = [l.split(":") for l in idx.strip().splitlines()]
    out: dict[str, str] = {}
    for n, parts in enumerate(lines):
        key = f"{parts[3]}:{parts[4]}"
        for name, want in FIELDS.items():
            if key == want and name not in out:
                start = int(parts[1])
                end = int(lines[n + 1][1]) - 1 if n + 1 < len(lines) else ""
                out[name] = f"{start}-{end}"
    missing = set(FIELDS) - set(out)
    if missing:
        raise RuntimeError(f"fields missing from idx: {missing}")
    return out


def decode(msg: bytes) -> tuple[np.ndarray, dict]:
    gid = eccodes.codes_new_from_message(msg)
    try:
        ni, nj = eccodes.codes_get(gid, "Ni"), eccodes.codes_get(gid, "Nj")
        vals = eccodes.codes_get_values(gid).reshape(nj, ni)
        geo = {k: eccodes.codes_get(gid, k) for k in (
            "Latin1InDegrees", "Latin2InDegrees", "LoVInDegrees", "LaDInDegrees",
            "latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees",
            "DxInMetres", "DyInMetres")}
    finally:
        eccodes.codes_release(gid)
    return vals, geo


def target_grid(geo: dict) -> tuple[np.ndarray, np.ndarray, int]:
    """Fractional (row, col) HRRR indices for every output pixel (Web Mercator rows)."""
    lon_0 = geo["LoVInDegrees"] - 360 if geo["LoVInDegrees"] > 180 else geo["LoVInDegrees"]
    p = Proj(proj="lcc", lat_1=geo["Latin1InDegrees"], lat_2=geo["Latin2InDegrees"],
             lat_0=geo["LaDInDegrees"], lon_0=lon_0, R=6371229)
    lon1 = geo["longitudeOfFirstGridPointInDegrees"]
    x0, y0 = p(lon1 - 360 if lon1 > 180 else lon1, geo["latitudeOfFirstGridPointInDegrees"])
    merc = lambda lat: math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))
    ymin, ymax = merc(SOUTH), merc(NORTH)
    height = round(WIDTH * (ymax - ymin) / math.radians(EAST - WEST))
    lons = np.linspace(WEST, EAST, WIDTH)
    ys = np.linspace(ymax, ymin, height)  # top row = north
    lats = np.degrees(2 * np.arctan(np.exp(ys)) - math.pi / 2)
    LON, LAT = np.meshgrid(lons, lats)
    x, y = p(LON, LAT)
    return (y - y0) / geo["DyInMetres"], (x - x0) / geo["DxInMetres"], height


def colorize(v: np.ndarray, scale, mask=None, alpha=235) -> np.ndarray:
    rgba = np.zeros(v.shape + (4,), np.uint8)
    for thr, hx in scale:
        sel = v >= thr if mask is None else (v >= thr) & mask
        rgba[sel] = (*hexrgb(hx), alpha)
    return rgba


def colorize_smooth(v: np.ndarray, stops, alpha_stops) -> np.ndarray:
    """Continuous colormap: linear interpolation between (value, hex) stops."""
    xs = [s[0] for s in stops]
    rgba = np.zeros(v.shape + (4,), np.uint8)
    for c in range(3):
        rgba[..., c] = np.interp(v, xs, [hexrgb(s[1])[c] for s in stops]).astype(np.uint8)
    a = np.interp(v, [s[0] for s in alpha_stops], [s[1] for s in alpha_stops])
    rgba[..., 3] = np.where(v >= alpha_stops[0][0], a, 0).astype(np.uint8)
    return rgba


def stp_fixed(cape, cin, srh1, shr6, t2, td2):
    """SPC fixed-layer significant tornado parameter (Thompson et al. 2003/2012)."""
    lcl = 125.0 * (t2 - td2)  # m, Espy approximation
    lcl_t = np.clip((2000 - lcl) / 1000, 0, 1)
    shr_t = np.where(shr6 < 12.5, 0, np.minimum(shr6, 30) / 20)
    cin_t = np.clip((200 + cin) / 150, 0, 1)
    return (cape / 1500) * lcl_t * (srh1 / 150) * shr_t * cin_t


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--hours", type=int, default=18)
    ap.add_argument("--run")
    a = ap.parse_args()
    run = dt.datetime.strptime(a.run, "%Y%m%d%H") if a.run else latest_run(a.hours)
    print(f"run {run:%Y-%m-%d %HZ}", flush=True)

    def download_once(fh: int) -> dict[str, bytes]:
        rngs = byte_ranges(fetch(url_for(run, fh, ".idx")).decode())
        return {name: fetch(url_for(run, fh), r) for name, r in rngs.items()}

    def download(fh: int) -> dict[str, bytes] | None:
        try:
            return with_retries(download_once, fh)
        except (urllib.error.URLError, OSError) as e:
            print(f"  f{fh:02d} download failed after retries: {e}", file=sys.stderr, flush=True)
            return None

    # Network in parallel; ecCodes is not thread-safe, so decode serially.
    with ThreadPoolExecutor(max_workers=6) as ex:
        raw = list(ex.map(download, range(1, a.hours + 1)))
    # Keep only the contiguous run of good hours; stop at the first failure.
    if None in raw:
        raw = raw[:raw.index(None)]
    min_hours = min(MIN_HOURS, a.hours)
    if len(raw) < min_hours:
        print(f"only {len(raw)} forecast hours downloaded (need {min_hours})", file=sys.stderr)
        return 1
    hours = []
    for msgs in raw:
        d = {}
        for name, msg in msgs.items():
            d[name], d["_geo"] = decode(msg)
        hours.append(d)

    rows, cols, height = target_grid(hours[0]["_geo"])
    lin = lambda f: map_coordinates(f, [rows, cols], order=1, mode="nearest")
    near = lambda f: map_coordinates(f, [rows, cols], order=0, mode="nearest")

    # Only the current run's frames: clear old product dirs so a longer previous run
    # can't leave stale fNN.png files behind.
    for p in ("refl", "stp", "uh"):
        shutil.rmtree(os.path.join(a.out, p), ignore_errors=True)
        os.makedirs(os.path.join(a.out, p))
    uh_max = None
    for fh, d in enumerate(hours, start=1):
        # Reflectivity + precip type
        refl = gaussian_filter(lin(d["refd"]), 1.0)
        snow, frzr, sleet = near(d["csnow"]) > 0, near(d["cfrzr"]) > 0, near(d["cicep"]) > 0
        rain = ~(snow | frzr | sleet)
        img = colorize(refl, RAIN, rain)
        for scale, m in ((SNOW, snow & ~frzr & ~sleet), (SLEET, sleet & ~frzr), (FRZR, frzr)):
            layer = colorize(refl, scale, m)
            img[layer[..., 3] > 0] = layer[layer[..., 3] > 0]
        Image.fromarray(img).save(os.path.join(a.out, "refl", f"f{fh:02d}.png"), optimize=True)
        # STP
        shr6 = np.hypot(d["ushr6"], d["vshr6"])
        stp = stp_fixed(d["cape"], d["cin"], d["srh1"], shr6, d["t2"], d["td2"])
        stp = gaussian_filter(lin(stp), 1.0)
        Image.fromarray(colorize(stp, STP, alpha=200)).save(
            os.path.join(a.out, "stp", f"f{fh:02d}.png"), optimize=True)
        # UH 2-5 km, running max through this hour
        uh_max = d["uh25"] if uh_max is None else np.maximum(uh_max, d["uh25"])
        Image.fromarray(colorize_smooth(lin(uh_max), UH_STOPS, UH_ALPHA)).save(
            os.path.join(a.out, "uh", f"f{fh:02d}.png"), optimize=True)
        print(f"  f{fh:02d} done", flush=True)

    meta = {"run": run.strftime("%Y-%m-%dT%H:00:00Z"), "hours": len(hours),
            "bounds": [[SOUTH, WEST], [NORTH, EAST]], "width": WIDTH, "height": height,
            "generated": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    # meta.json last and atomically, so a reader never sees meta pointing at missing frames.
    tmp = os.path.join(a.out, "meta.json.tmp")
    with open(tmp, "w") as f:
        json.dump(meta, f)
    os.replace(tmp, os.path.join(a.out, "meta.json"))
    print("ok", meta)
    return 0


if __name__ == "__main__":
    sys.exit(main())
