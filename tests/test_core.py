import io
import math
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import pytest

import adsb_coverage as ac

RE = 4 / 3 * ac.core.R_EARTH
FT = ac.core.FT


class Flat:
    def sample(self, lat, lon):
        return np.zeros(lat.shape)


class Ridge:
    """Crête de 2000 m à 20 km à l'est de (46.9, 7.4)."""

    def sample(self, lat, lon):
        dx = (lon - 7.4) * 111320 * math.cos(math.radians(46.9))
        dy = (lat - 46.9) * 111200
        return 2000 * np.exp(-((dx - 20000) ** 2) / (2 * 1500**2)) * np.exp(-(dy**2) / (2 * 8000**2))


def test_lv95_berne():
    e, n = ac.wgs84_to_lv95(46.9510827861, 7.4386324167)  # ancien observatoire de Berne
    assert abs(e - 2600000) < 2 and abs(n - 1200000) < 2


@pytest.mark.parametrize("alt_m", [305.0, 914.0, 9144.0])
def test_flat_amsl_matches_radio_horizon(alt_m):
    _, r, _ = ac.coverage_ranges(ac.Terrain([Flat()]), 46.9, 7.4, 10.0, [alt_m], 450e3, 100.0, 10.0, 4 / 3)
    assert abs(r[0].mean() - (math.sqrt(2 * RE * 10) + math.sqrt(2 * RE * alt_m))) < 250


@pytest.mark.parametrize("h_ft", [100, 300, 1000])
def test_flat_agl_matches_radio_horizon(h_ft):
    _, r, _ = ac.coverage_ranges(ac.Terrain([Flat()]), 46.9, 7.4, 10.0, [h_ft * FT], 150e3, 100.0, 10.0, 4 / 3, agl=True)
    assert abs(r[0].mean() - (math.sqrt(2 * RE * 10) + math.sqrt(2 * RE * h_ft * FT))) < 250


def test_ridge_shadow():
    az, r, _ = ac.coverage_ranges(ac.Terrain([Ridge()]), 46.9, 7.4, 10.0, [3000.0], 400e3, 100.0, 10.0, 4 / 3)
    assert r[0][list(az).index(90)] < 40e3 < 150e3 < r[0][list(az).index(270)]


def test_swissalti3d_xyz_decimation_and_sampling():
    rows = ["X Y Z"] + [f"{2600000 + i * 2 + 1}.0 {1200999 - j * 2}.0 {500 + i * 0.2:.2f}"
                        for j in range(500) for i in range(500)]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("tile.xyz", "\n".join(rows))
    sw = ac.SwissAlti3D.__new__(ac.SwissAlti3D)
    sw.f, sw.res, sw.n, sw.tiles = 25, 50.0, 20, {}
    sw.tiles[2600 * 100000 + 1200] = sw._decimate(buf.getvalue())
    z = sw.sample(np.array([46.9555]), np.array([7.4452]))  # ~ E 2600500, N 1200500
    assert 500 <= z[0] <= 600


def test_cli_kml_valid(tmp_path):
    out = tmp_path / "c.kml"
    ac.main(["46.9", "7.4", "--flat", "--agl", "--az-step", "10", "-o", str(out)])
    root = ET.parse(out).getroot()
    assert len(root.findall(".//{http://www.opengis.net/kml/2.2}Polygon")) > 0
