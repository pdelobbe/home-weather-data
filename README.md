# home-weather-data

Hourly HRRR map overlays for the Home Weather app (central Oklahoma): 1 km reflectivity colored by precipitation type, significant tornado parameter, and 2–5 km updraft helicity (run max).

- `.github/workflows/hrrr.yml` renders the newest complete HRRR run from NOAA's open data on AWS (`noaa-hrrr-bdp-pds`) twice an hour.
- Output is force-pushed to the `wx-data` branch as one commit (`meta.json`, `refl/fNN.png`, `stp/fNN.png`, `uh/fNN.png`), so history never grows.
- `hrrr_render.py` is a copy of `pipeline/hrrr_render.py` from the app repo — edit it there and copy it here.

Data: NOAA HRRR (public domain).
