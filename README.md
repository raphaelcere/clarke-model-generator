# adsb-coverage

**Theoretical** reception coverage of a fixed ADS-B (1090 MHz) / ADS-L ground receiver: line of sight (LOS)
with an **effective Earth radius** (k = 4/3 by default) and **terrain masking**. Input: latitude and longitude.
Output: a **KML** file (Google Earth, QGIS).

## At a glance (for a general audience)

An ADS-B / ADS-L receiver picks up the messages that aircraft broadcast continuously (position, altitude,
identification). These radio waves travel **in a straight line**: if a mountain lies between the aircraft and
the antenna, the signal is blocked, and the curvature of the Earth also hides aircraft that are low and far away.
The tool computes **how far** a given receiver can "see" aircraft, depending on their altitude.

![Principle: an aircraft is received if it is in direct line of sight of the antenna](docs/principle.svg)

The result is a map of the reception zone, one zone per flight level: the higher an aircraft flies, the farther
it can be seen, and the terrain hollows out "holes" behind mountains, especially for low flights.

![Reception zone around the receiver, seen from above](docs/coverage-zone.svg)

## Installation

```
pip install -e .            # numpy + requests
pip install -e ".[geotiff]" # + rasterio (far terrain from Copernicus, local terrain models)
```

Python 3.9 to 3.13 (tested by the continuous integration on every push).

## Usage

```
# altitudes AMSL (default: 1000 ... 30000 ft)
adsb-coverage 46.9480 7.4474 --ant-h 15 -o bern.kml

# low altitudes, heights above the local ground
adsb-coverage 46.9480 7.4474 --agl --alt-ft 30,100,300,1000 --near-res-m 10 -o low.kml
```

A summary of the ranges is printed on the console, and the KML has one layer per level (with `--agl`, the visible
sectors are drawn one by one). `adsb-coverage --help` lists all the options.

## Online calculation (GitHub Actions)

The repository contains a workflow that computes a coverage without installing anything, from GitHub's servers
(direct access to swisstopo and Copernicus, no corporate proxy):

1. Open the **Actions** tab, select the workflow **"ADS-B / ADS-L coverage"** and click **Run workflow**.
2. Enter the latitude, longitude, mode (`agl` or `amsl`), antenna height and, if needed, the levels in ft.
3. When the run has finished, download the KML from the **Artifacts** section (kept for 30 days); the table of
   ranges appears in the run summary.

The inputs are validated (numbers only, bounded ranges) before anything is computed. Write access to the
repository is needed to start the workflow. On a public repository, the parameters entered are visible in the
run history.

## Example

[`examples/bern-agl-low-altitudes.kml`](examples/bern-agl-low-altitudes.kml): low-altitude coverage
(`--agl` mode) computed by the GitHub Actions workflow for a receiver at 46.9744 N, 7.4771 E
(ground at 521 m, antenna 10 m above ground, k = 4/3), for aircraft 15, 30, 60, 90, 120 and 150 ft above ground.

| Height above ground | Max continuous range | Max visible range |
|---|---|---|
| 150 ft | 6 km | 82 km |
| 120 ft | 5 km | 81 km |
| 90 ft | 4 km | 81 km |
| 60 ft | 4 km | 81 km |
| 30 ft | 2 km | 81 km |
| 15 ft | 2 km | 81 km |

*Continuous range*: distance to the first obstruction by the terrain. *Visible range*: farthest visible point, even
behind a ridge (at low altitude an aircraft can reappear above a distant, higher relief).

**Opening it in Google Earth**

1. **Google Earth Pro** (desktop): [download the file](https://github.com/raphaelcere/clarke-model-generator/raw/main/examples/bern-agl-low-altitudes.kml)
   (right-click, "Save link as..."), then double-click the `.kml`, or *File, Open*.
2. **Google Earth web** ([earth.google.com/web](https://earth.google.com/web/)): *Projects*, *Open*, *Import KML file from computer*,
   then pick the downloaded file. The file is large (about 1 MB), so the import can take a few seconds.
3. **QGIS**: drag and drop the file into the window.

## Architecture (UML)

[Mermaid](https://mermaid.js.org/) diagrams, rendered directly by GitHub. `TerrainSource` is a convention
(duck typing: every terrain source exposes `sample(lat, lon)`), not a class in the code.

### Classes and modules

```mermaid
classDiagram
    direction LR
    class TerrainSource {
        <<interface>>
        +sample(lat, lon) elevations
    }
    class SwissAlti3D {
        +res : float
        +tiles : dict
        -_load(lat, lon, radius)
        -_download(...)
        -_list_items(bbox)
        -_fetch(tile)
        +sample(lat, lon)
    }
    class GeoTiffMosaic {
        +chunks : list
        +sample(lat, lon)
    }
    class Terrain {
        +sources : list
        +unknown : float
        +sample(lat, lon)
    }
    class Propagation {
        <<module functions>>
        +compute_coverage(...) visible
        +continuous_range(visible, d)
        +max_visible_range(visible, d)
        +visible_runs(visible, d, step)
    }
    class KmlWriter {
        <<module functions>>
        +write_kml(...)
    }
    class CLI {
        <<main>>
        +main(argv)
    }
    class Geodesy {
        <<module functions>>
        +wgs84_to_lv95(lat, lon)
        +destination(lat, lon, az, dist)
    }
    TerrainSource <|.. SwissAlti3D
    TerrainSource <|.. GeoTiffMosaic
    TerrainSource <|.. Terrain
    Terrain o-- "1..*" TerrainSource : tries in order
    CLI --> Terrain : builds
    CLI --> Propagation : calls
    CLI --> KmlWriter : calls
    Propagation --> Terrain : samples the terrain
    Propagation ..> Geodesy
    SwissAlti3D ..> Geodesy
    KmlWriter ..> Geodesy
```

### Sequence of a calculation

```mermaid
sequenceDiagram
    actor U as User
    participant C as CLI (main)
    participant T as Terrain
    participant S as SwissAlti3D
    participant G as GeoTiffMosaic
    participant P as compute_coverage
    participant K as write_kml
    U->>C: adsb-coverage lat lon [--agl] -o file.kml
    C->>S: load the 2 m tiles (near zone)
    S-->>C: decimated terrain (local cache)
    C->>G: read Copernicus 30 m or a local model (far zone)
    C->>T: chain the sources
    C->>P: antenna, levels, maximum range, k
    loop for each azimuth and distance step
        P->>T: ground elevation?
        T-->>P: terrain (swissALTI3D, else Copernicus)
    end
    P-->>C: visibility mask per level
    C->>K: polygons per level
    K-->>U: KML file + table of ranges
```

## Model

For each azimuth, the terrain is sampled along a profile. In a "flat Earth" frame (apparent height
h - d^2/(2 k R)), an aircraft is visible if the slope from the antenna to the aircraft is at least the steepest slope
to any terrain met so far. There is no Fresnel zone, link budget or man-made obstacle (buildings, vegetation):
this is purely geometric visibility.

## Terrain data

| Source | Use | Notes |
|---|---|---|
| [swissALTI3D](https://www.swisstopo.admin.ch/en/height-model-swissalti3d) (swisstopo, 2 m) | near zone (`--near-km`, default 10 km) | downloaded automatically through the STAC API, cached in `~/.cache/adsb_coverage`, decimated to `--near-res-m` (default 50 m). Bare ground. |
| [Copernicus GLO-30](https://registry.opendata.aws/copernicus-dem/) (30 m) | rest of the domain | low-resolution read (`--far-res-m`, default 300 m), requires `rasterio` |
| Local EPSG:4326 GeoTIFF | `--dem file.tif` | |

Corporate network: `--proxy http://host:port`, `--ca-bundle ca.pem`, or offline mode with `--tile-dir`
(swissALTI3D `.xyz.zip` tiles downloaded by hand).

## Known limitations and to do

- Copernicus tile download and GeoTIFF reading still need validating in the field (the continuous integration only
  tests a synthetic GeoTIFF).
- Migration of the STAC API from v0.9 (deprecated, active at least until the end of 2026) to v1.
- `--alt-m` option, automatic detection of the system proxy.
- Comparison with the range observed in real logs (per azimuth), calibration of k.
- Man-made obstacles (swissSURFACE3D), Fresnel zone, link budget.

## Development

```
pip install -e ".[dev]"
pytest
```

The tests are also run by GitHub Actions (workflow "Tests") on Python 3.9 to 3.13 at every `push` to `main`.

Licence: to be defined.
