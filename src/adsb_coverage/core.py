"""
adsb_coverage.py - Couverture de réception théorique d'un récepteur fixe ADS-B / ADS-L (au sol).

Modèle
------
Visibilité radioélectrique directe (LOS) avec rayon terrestre effectif k·R (k = 4/3 par défaut,
modèle d'atmosphère standard « de Clarke »), obstruction par le relief le long de rayons radiaux.
Pour chaque azimut et chaque altitude d'aéronef (AMSL), on calcule la portée maximale à laquelle
un aéronef reste visible depuis l'antenne.

Topographie
-----------
- swissALTI3D (swisstopo, 2 m) téléchargé automatiquement, pour la zone proche (--near-km).
- Copernicus GLO-30 (30 m, AWS open data) pour le reste (hors Suisse / champ lointain), lu en
  basse résolution via les overviews du COG. Nécessite rasterio.
- ou un MNT GeoTIFF local en EPSG:4326 (--dem).

Sortie : KML (un polygone par altitude) + résumé des portées en console.

Exemples
--------
    adsb-coverage 46.9480 7.4474 --ant-h 15 --alt-ft 1000,3000,10000,30000 -o bern.kml
    adsb-coverage 46.9480 7.4474 --agl -o bas.kml        # basses altitudes, hauteurs sol

Dépendances : numpy, requests ; rasterio (champ lointain Copernicus / MNT local).
"""
import argparse
import io
import math
import sys
import warnings
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np

R_EARTH = 6371008.8  # m, rayon moyen
FT = 0.3048
STAC_ITEMS = "https://data.geo.admin.ch/api/stac/v0.9/collections/ch.swisstopo.swissalti3d/items"
COP_URL = "https://copernicus-dem-30m.s3.amazonaws.com/{n}/{n}.tif"


# ----------------------------------------------------------------------------- réseau
HTTP = {"proxies": None, "verify": True, "timeout": 30}


def http_get(url, **kw):
    """requests.get avec proxy / CA configurables (--proxy, --ca-bundle ; sinon variables HTTPS_PROXY, etc.)."""
    import requests

    kw.setdefault("timeout", HTTP["timeout"])
    return requests.get(url, proxies=HTTP["proxies"], verify=HTTP["verify"], **kw)


def network_help(exc):
    return (f"Connexion impossible : {exc}\n"
            "Pistes : (1) proxy -> --proxy http://utilisateur:mdp@hote:port ou variable HTTPS_PROXY ; "
            "(2) inspection SSL -> --ca-bundle chemin/vers/ca.pem ; "
            "(3) mode hors ligne -> télécharger les tuiles swissALTI3D (2 m, .xyz.zip) depuis un poste avec accès "
            "et les placer dans un dossier passé via --tile-dir ; (4) --near-km 0 --dem mnt.tif.")


# ----------------------------------------------------------------------------- géodésie
def wgs84_to_lv95(lat, lon):
    """WGS84 -> MN95/LV95 (formules approchées swisstopo, ~1 m en Suisse)."""
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    p = (lat * 3600 - 169028.66) / 10000
    l = (lon * 3600 - 26782.5) / 10000
    e = 2600072.37 + 211455.93 * l - 10938.51 * l * p - 0.36 * l * p**2 - 44.54 * l**3
    n = 1200147.07 + 308807.95 * p + 3745.25 * l**2 + 76.63 * p**2 - 194.56 * l**2 * p + 119.79 * p**3
    return e, n


def destination(lat, lon, az_deg, dist_m):
    """Point atteint depuis (lat, lon) à l'azimut et la distance donnés (sphère). Broadcast numpy."""
    p1, l1 = math.radians(lat), math.radians(lon)
    th = np.radians(az_deg)
    d = dist_m / R_EARTH
    p2 = np.arcsin(np.sin(p1) * np.cos(d) + np.cos(p1) * np.sin(d) * np.cos(th))
    l2 = l1 + np.arctan2(np.sin(th) * np.sin(d) * np.cos(p1), np.cos(d) - np.sin(p1) * np.sin(p2))
    return np.degrees(p2), (np.degrees(l2) + 540.0) % 360.0 - 180.0


# ----------------------------------------------------------------------------- MNT : swissALTI3D
class SwissAlti3D:
    """Tuiles 1 km swissALTI3D (2 m), décimées à `res_m` et mises en cache en .npy."""

    def __init__(self, lat, lon, radius_m, res_m=50, cache_dir=None, workers=8, tile_dir=None):
        valid = [f for f in (1, 2, 4, 5, 10, 20, 25, 50) if 500 % f == 0]
        self.f = min(valid, key=lambda f: abs(f * 2 - res_m))
        self.res = self.f * 2.0
        self.n = 500 // self.f
        self.cache = Path(cache_dir or Path.home() / ".cache" / "adsb_coverage" / "swissalti3d")
        self.cache.mkdir(parents=True, exist_ok=True)
        self.workers = workers
        self.tile_dir = Path(tile_dir) if tile_dir else None
        self.tiles = {}
        self._load(lat, lon, radius_m)

    def _wanted(self, lat, lon, r):
        e0, n0 = wgs84_to_lv95(lat, lon)
        out = []
        for tx in range(int((e0 - r) // 1000), int((e0 + r) // 1000) + 1):
            for ty in range(int((n0 - r) // 1000), int((n0 + r) // 1000) + 1):
                cx, cy = tx * 1000 + 500, ty * 1000 + 500
                if math.hypot(cx - e0, cy - n0) <= r + 710:
                    out.append((tx, ty))  # id swisstopo : E/1000 - N/1000, ex. 2485-1109
        return out

    def _path(self, tx, ty):
        return self.cache / f"{tx}-{ty}_{int(self.res)}m.npy"

    def _list_items(self, bbox):
        url, params, out = STAC_ITEMS, {"bbox": ",".join(f"{v:.6f}" for v in bbox), "limit": 100}, {}
        while url:
            r = http_get(url, params=params, timeout=60)
            r.raise_for_status()
            j = r.json()
            for it in j.get("features", []):
                _, year, key = it["id"].split("_")
                href = next((a["href"] for k, a in it["assets"].items()
                             if "_2_2056" in k and k.endswith(".xyz.zip")), None)
                if href and (key not in out or year > out[key][0]):
                    out[key] = (year, href)
            url = next((l["href"] for l in j.get("links", []) if l.get("rel") == "next"), None)
            params = None
        return out

    def _decimate(self, zbytes):
        with zipfile.ZipFile(io.BytesIO(zbytes)) as z:
            name = next(n for n in z.namelist() if n.lower().endswith(".xyz"))
            text = z.read(name).decode("ascii", "ignore")
        lines = text.splitlines()
        if lines and not lines[0].lstrip()[:1].isdigit():
            lines = lines[1:]  # en-tête « X Y Z »
        a = np.array(" ".join(lines).split(), dtype=float).reshape(-1, 3)
        ix = np.rint((a[:, 0] - a[:, 0].min()) / 2.0).astype(int)
        iy = np.rint((a[:, 1].max() - a[:, 1]) / 2.0).astype(int)  # ligne 0 = nord
        g = np.full((500, 500), np.nan)
        g[np.clip(iy, 0, 499), np.clip(ix, 0, 499)] = a[:, 2]
        g[g < -500] = np.nan
        f = self.f
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmean(g.reshape(500 // f, f, 500 // f, f), axis=(1, 3)).astype(np.float32)

    def _fetch(self, tile):
        (tx, ty), (year, href) = tile
        for attempt in range(3):
            try:
                r = http_get(href, timeout=120)
                r.raise_for_status()
                arr = self._decimate(r.content)
                np.save(self._path(tx, ty), arr)
                return (tx, ty), arr
            except Exception as exc:  # noqa: BLE001
                err = exc
        print(f"  ! tuile {tx}-{ty} ({year}) échec : {err}", file=sys.stderr)
        return (tx, ty), None

    def _load(self, lat, lon, r):
        want = self._wanted(lat, lon, r)
        missing = []
        for tx, ty in want:
            p = self._path(tx, ty)
            if p.exists():
                self.tiles[tx * 100000 + ty] = np.load(p)
            else:
                missing.append((tx, ty))
        if self.tile_dir and missing:  # tuiles .xyz.zip fournies manuellement
            still = []
            for tx, ty in missing:
                f = next(iter(sorted(self.tile_dir.glob(f"*{tx}-{ty}_2_2056*.xyz.zip"))), None)
                if f is None:
                    still.append((tx, ty))
                    continue
                arr = self._decimate(f.read_bytes())
                np.save(self._path(tx, ty), arr)
                self.tiles[tx * 100000 + ty] = arr
            missing = still
            if missing:
                print(f"swissALTI3D : {len(missing)} tuile(s) absentes de --tile-dir", file=sys.stderr)
                return
        if not missing:
            return
        print(f"swissALTI3D : {len(want)} tuiles, {len(missing)} à télécharger…", file=sys.stderr)
        dlat, dlon = r / 111_200, r / (111_200 * math.cos(math.radians(lat)))
        try:
            items = self._list_items((lon - dlon, lat - dlat, lon + dlon, lat + dlat))
        except Exception as exc:  # noqa: BLE001
            raise SystemExit(network_help(exc)) from None
        todo = [((tx, ty), items[f"{tx}-{ty}"]) for tx, ty in missing if f"{tx}-{ty}" in items]
        with ThreadPoolExecutor(self.workers) as ex:
            for (tx, ty), arr in ex.map(self._fetch, todo):
                if arr is not None:
                    self.tiles[tx * 100000 + ty] = arr
        print(f"swissALTI3D : {len(self.tiles)}/{len(want)} tuiles disponibles", file=sys.stderr)

    def sample(self, lat, lon):
        e, n = wgs84_to_lv95(lat, lon)
        out = np.full(e.shape, np.nan)
        tx, ty = np.floor(e / 1000).astype(np.int64), np.floor(n / 1000).astype(np.int64)
        keys, inv = np.unique(tx * 100000 + ty, return_inverse=True)
        inv = inv.reshape(e.shape)
        for i, k in enumerate(keys):
            t = self.tiles.get(int(k))
            if t is None:
                continue
            m = inv == i
            col = np.clip(((e[m] - (k // 100000) * 1000) / self.res).astype(int), 0, self.n - 1)
            row = np.clip((((k % 100000) * 1000 + 1000 - n[m]) / self.res).astype(int), 0, self.n - 1)
            out[m] = t[row, col]
        return out


# ----------------------------------------------------------------------------- MNT : GeoTIFF (WGS84)
class GeoTiffMosaic:
    """Mosaïque de GeoTIFF en EPSG:4326 lus en basse résolution (overviews des COG)."""

    def __init__(self, paths, bbox, res_m=300):
        import rasterio
        from rasterio.enums import Resampling
        from rasterio.windows import Window

        lon0, lat0, lon1, lat1 = bbox
        self.chunks = []
        env = rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif")
        with env:
            for p in paths:
                try:
                    ds = rasterio.open(p)
                except Exception:  # noqa: BLE001  (tuile absente : mer, hors couverture)
                    continue
                with ds:
                    if not (ds.crs and ds.crs.is_geographic):
                        raise SystemExit(f"{p} : le MNT doit être en coordonnées géographiques (EPSG:4326).")
                    b, t = ds.bounds, ds.transform
                    w, e, s, n = max(lon0, b.left), min(lon1, b.right), max(lat0, b.bottom), min(lat1, b.top)
                    if w >= e or s >= n:
                        continue
                    c0, r0 = (~t) * (w, n)
                    c1, r1 = (~t) * (e, s)
                    c0, r0 = max(int(math.floor(c0)), 0), max(int(math.floor(r0)), 0)
                    c1, r1 = min(int(math.ceil(c1)), ds.width), min(int(math.ceil(r1)), ds.height)
                    step_deg = res_m / 111_320
                    ow = max(1, int(round((c1 - c0) * abs(t.a) / step_deg)))
                    oh = max(1, int(round((r1 - r0) * abs(t.e) / step_deg)))
                    win = Window(c0, r0, c1 - c0, r1 - r0)
                    try:
                        arr = ds.read(1, window=win, out_shape=(oh, ow), resampling=Resampling.max)
                    except Exception:  # noqa: BLE001
                        arr = ds.read(1, window=win, out_shape=(oh, ow), resampling=Resampling.average)
                    arr = arr.astype(np.float32)
                    if ds.nodata is not None:
                        arr[arr == ds.nodata] = np.nan
                    self.chunks.append((t.c + c0 * t.a, t.f + r0 * t.e, (c1 - c0) * t.a / ow,
                                        (r1 - r0) * abs(t.e) / oh, arr))
        print(f"MNT GeoTIFF : {len(self.chunks)} tuile(s) lue(s)", file=sys.stderr)

    def sample(self, lat, lon):
        out = np.full(lat.shape, np.nan)
        for west, north, dx, dy, arr in self.chunks:
            col = np.floor((lon - west) / dx).astype(int)
            row = np.floor((north - lat) / dy).astype(int)
            m = (col >= 0) & (col < arr.shape[1]) & (row >= 0) & (row < arr.shape[0]) & np.isnan(out)
            out[m] = arr[row[m], col[m]]
        return out


def copernicus_urls(bbox):
    lon0, lat0, lon1, lat1 = bbox
    urls = []
    for la in range(math.floor(lat0), math.floor(lat1) + 1):
        for lo in range(math.floor(lon0), math.floor(lon1) + 1):
            n = (f"Copernicus_DSM_COG_10_{'N' if la >= 0 else 'S'}{abs(la):02d}_00_"
                 f"{'E' if lo >= 0 else 'W'}{abs(lo):03d}_00_DEM")
            urls.append(COP_URL.format(n=n))
    return urls


class Terrain:
    """Chaîne de MNT : le premier qui renvoie une valeur gagne ; sinon 0 m (mer / inconnu)."""

    def __init__(self, sources):
        self.sources = [s for s in sources if s is not None]

    def sample(self, lat, lon):
        lat, lon = np.asarray(lat), np.asarray(lon)
        z = np.full(lat.shape, np.nan)
        for s in self.sources:
            m = np.isnan(z)
            if not m.any():
                break
            z[m] = s.sample(lat[m], lon[m])
        self.unknown = float(np.isnan(z).mean())
        return np.nan_to_num(z, nan=0.0)


# ----------------------------------------------------------------------------- modèle de propagation
def compute_coverage(terrain, lat, lon, ant_h, heights_m, rmax, ds, az_step, k, ground=None, agl=False):
    """
    Visibilité de l'antenne le long de rayons radiaux -> masque booléen vis[alt, azimut, distance].

    Repère « terre aplatie » : un point à la distance d et à la hauteur h a la hauteur apparente
    h - d²/(2·k·R). Un aéronef est visible si la pente antenne->aéronef est >= à la pente maximale
    du relief rencontré jusque-là (relief au point lui-même inclus).
    - agl=False : hauteurs = altitudes AMSL constantes ; la visibilité est monotone avec la distance.
    - agl=True  : hauteurs = hauteur au-dessus du sol local (l'aéronef suit le relief) ; la visibilité
      n'est plus monotone (un aéronef peut réapparaître au-delà d'une crête, dans une vallée plus haute).
    """
    re_eff = k * R_EARTH
    az = np.arange(0, 360, az_step, dtype=float)
    d = np.arange(ds, rmax + ds, ds, dtype=float)
    plat, plon = destination(lat, lon, az[:, None], d[None, :])
    z0 = float(terrain.sample(np.array([lat]), np.array([lon]))[0]) if ground is None else ground
    z = terrain.sample(plat, plon)
    h0 = z0 + ant_h
    bulge = d**2 / (2 * re_eff)
    m_terr = np.maximum.accumulate((z - bulge - h0) / d[None, :], axis=1)
    vis = np.empty((len(heights_m), len(az), len(d)), dtype=bool)
    for i, h in enumerate(heights_m):
        top = (z + h) if agl else h  # altitude AMSL de l'aéronef : (az, d) ou scalaire
        slope_ac = (top - bulge[None, :] - h0) / d[None, :]
        vis[i] = slope_ac >= m_terr
    return az, d, vis, z0


def continuous_range(vis, d):
    """Distance jusqu'au premier masquage (couverture continue depuis le récepteur), par azimut."""
    blocked = ~vis
    first = np.where(blocked.any(axis=-1), blocked.argmax(axis=-1), vis.shape[-1])
    return np.where(first > 0, d[np.clip(first - 1, 0, len(d) - 1)], 0.0)


def max_visible_range(vis, d):
    """Distance du point visible le plus lointain, par azimut (0 si aucun)."""
    last = vis.shape[-1] - 1 - vis[..., ::-1].argmax(axis=-1)
    return np.where(vis.any(axis=-1), d[last], 0.0)


def coverage_ranges(terrain, lat, lon, ant_h, alts_m, rmax, ds, az_step, k, ground=None, agl=False):
    """Compatibilité : (azimuts, portée continue [alt, azimut] en m, altitude du sol au récepteur)."""
    az, d, vis, z0 = compute_coverage(terrain, lat, lon, ant_h, alts_m, rmax, ds, az_step, k, ground, agl)
    return az, np.stack([continuous_range(v, d) for v in vis]), z0


def visible_runs(vis_alt, d, ds, min_len=2):
    """Segments visibles contigus par azimut : (indice azimut, d_debut, d_fin) en m, longueur >= min_len pas."""
    pad = np.zeros((vis_alt.shape[0], 1), dtype=int)
    diff = np.diff(np.concatenate([pad, vis_alt.astype(int), pad], axis=1), axis=1)
    (ia, st), (_, en) = np.nonzero(diff == 1), np.nonzero(diff == -1)
    keep = (en - st) >= min_len
    ia, st, en = ia[keep], st[keep], en[keep]
    return ia, d[st] - ds / 2, d[en - 1] + ds / 2


# ----------------------------------------------------------------------------- KML
def _colors(n):
    for i in range(n):
        t = i / max(n - 1, 1)
        yield int(255 * t), int(200 * (1 - abs(2 * t - 1)) + 55), int(255 * (1 - t))


def write_kml(path, lat, lon, z0, ant_h, az, d, vis, heights_ft, k, name, agl=False, az_step=1.0, ds=250.0):
    n = len(heights_ft)
    ref = "AGL (hauteur sol)" if agl else "AMSL"
    cont = np.stack([continuous_range(v, d) for v in vis])
    vmax = np.stack([max_visible_range(v, d) for v in vis])
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<kml xmlns="http://www.opengis.net/kml/2.2"><Document>',
           f"<name>{escape(name)}</name>",
           f"<description>Couverture LOS théorique, k={k:.4f}, antenne {ant_h:g} m AGL "
           f"(sol {z0:.0f} m). Niveaux en {ref}.</description>"]
    for i, (r, g, b) in enumerate(_colors(n)):
        out.append(f'<Style id="s{i}"><LineStyle><color>ff{b:02x}{g:02x}{r:02x}</color><width>2</width></LineStyle>'
                   f'<PolyStyle><color>{"55" if agl else "33"}{b:02x}{g:02x}{r:02x}</color>'
                   f'<outline>{0 if agl else 1}</outline></PolyStyle></Style>')
    out.append('<Style id="rx"><IconStyle><scale>1.2</scale></IconStyle></Style>')
    out.append(f'<Placemark><name>Récepteur</name><styleUrl>#rx</styleUrl><description>Sol {z0:.0f} m, '
               f'antenne +{ant_h:g} m</description><Point><coordinates>{lon:.6f},{lat:.6f},0</coordinates></Point></Placemark>')
    out.append("<Folder><name>Couverture par niveau</name>")
    for i in reversed(range(n)):  # niveaux hauts dessous, bas dessus
        frac = vis[i].mean()
        label = (f"{heights_ft[i]:g} ft {ref} - continue max {cont[i].max() / 1000:.0f} km, "
                 f"visible max {vmax[i].max() / 1000:.0f} km")
        if agl:
            ia, d0, d1 = visible_runs(vis[i], d, ds)
            h = az_step / 2
            polys = []
            for a_i, r0, r1 in zip(ia, d0, d1):
                a0, a1 = az[a_i] - h, az[a_i] + h
                la = [destination(lat, lon, x, y) for x, y in ((a0, r0), (a1, r0), (a1, r1), (a0, r1))]
                ring = [(float(p[1]), float(p[0])) for p in la]
                ring.append(ring[0])
                polys.append("<Polygon><altitudeMode>clampToGround</altitudeMode><outerBoundaryIs><LinearRing><coordinates>"
                             + " ".join(f"{x:.5f},{y:.5f},0" for x, y in ring) + "</coordinates></LinearRing></outerBoundaryIs></Polygon>")
            if not polys:
                continue
            out.append(f'<Placemark><name>{escape(label)}</name><description>{frac:.1%} des points échantillonnés visibles'
                       f'</description><styleUrl>#s{i}</styleUrl><MultiGeometry>{"".join(polys)}</MultiGeometry></Placemark>')
        else:
            rr = cont[i]
            pl, pn = destination(lat, lon, az, rr)
            pl, pn = np.where(rr > 0, pl, lat), np.where(rr > 0, pn, lon)
            ring = list(zip(pn, pl)) + [(pn[0], pl[0])]
            coords = " ".join(f"{x:.6f},{y:.6f},0" for x, y in ring)
            out.append(f'<Placemark><name>{escape(label)}</name><styleUrl>#s{i}</styleUrl><Polygon><tessellate>1</tessellate>'
                       f'<altitudeMode>clampToGround</altitudeMode><outerBoundaryIs><LinearRing><coordinates>{coords}'
                       f'</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>')
    out.append("</Folder></Document></kml>")
    Path(path).write_text("\n".join(out), encoding="utf-8")
    return cont, vmax


# ----------------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lat", type=float, help="latitude WGS84 (degrés)")
    ap.add_argument("lon", type=float, help="longitude WGS84 (degrés)")
    ap.add_argument("-o", "--output", default="coverage.kml")
    ap.add_argument("--ant-h", type=float, default=10.0, help="hauteur de l'antenne au-dessus du sol (m)")
    ap.add_argument("--ground-m", type=float, default=None, help="altitude du sol au récepteur (sinon lue dans le MNT)")
    ap.add_argument("--agl", action="store_true",
                    help="niveaux en hauteur au-dessus du sol local (basses altitudes) au lieu d'altitudes AMSL")
    ap.add_argument("--alt-ft", default=None,
                    help="niveaux en ft, séparés par des virgules (défaut AMSL : 1000..30000 ; défaut AGL : 100..3000)")
    ap.add_argument("--k", type=float, default=4 / 3, help="facteur de rayon terrestre effectif (défaut 4/3)")
    ap.add_argument("--rmax-km", type=float, default=None, help="portée max calculée (défaut 400 AMSL, 150 AGL)")
    ap.add_argument("--step-m", type=float, default=None, help="pas le long des rayons (défaut 250 m AMSL, 100 m AGL)")
    ap.add_argument("--az-step", type=float, default=1.0, help="pas en azimut (degrés)")
    ap.add_argument("--near-km", type=float, default=10.0, help="rayon de la zone swissALTI3D (0 = désactivé)")
    ap.add_argument("--near-res-m", type=float, default=50.0, help="résolution après décimation de swissALTI3D (20 = plus fin)")
    ap.add_argument("--far", choices=["copernicus", "none"], default="copernicus", help="MNT hors zone proche")
    ap.add_argument("--far-res-m", type=float, default=300.0)
    ap.add_argument("--dem", nargs="*", default=None, help="GeoTIFF(s) EPSG:4326 locaux (remplacent Copernicus)")
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--tile-dir", default=None, help="dossier de tuiles swissALTI3D .xyz.zip téléchargées à la main")
    ap.add_argument("--proxy", default=None, help="proxy HTTP(S), ex. http://hote:8080")
    ap.add_argument("--ca-bundle", default=None, help="certificat(s) CA PEM (proxy avec inspection SSL)")
    ap.add_argument("--timeout", type=float, default=30.0, help="délai de connexion (s)")
    ap.add_argument("--flat", action="store_true", help="terrain plat à 0 m (horizon radio théorique, pour test)")
    a = ap.parse_args(argv)

    if a.proxy:
        HTTP["proxies"] = {"http": a.proxy, "https": a.proxy}
    if a.ca_bundle:
        HTTP["verify"] = a.ca_bundle
    HTTP["timeout"] = a.timeout
    alt_txt = a.alt_ft or ("100,300,500,1000,2000,3000" if a.agl else "1000,3000,6000,10000,20000,30000")
    heights_ft = [float(x) for x in alt_txt.split(",")]
    rmax = (a.rmax_km or (150.0 if a.agl else 400.0)) * 1000
    step = a.step_m or (100.0 if a.agl else 250.0)
    dlat = rmax / 111_200
    dlon = rmax / (111_200 * math.cos(math.radians(a.lat)))
    bbox = (a.lon - dlon, a.lat - dlat, a.lon + dlon, a.lat + dlat)

    sources = []
    if not a.flat:
        if a.near_km > 0:
            sources.append(SwissAlti3D(a.lat, a.lon, a.near_km * 1000, a.near_res_m, a.cache_dir, tile_dir=a.tile_dir))
        if a.dem:
            sources.append(GeoTiffMosaic(a.dem, bbox, a.far_res_m))
        elif a.far == "copernicus":
            sources.append(GeoTiffMosaic(copernicus_urls(bbox), bbox, a.far_res_m))
    terrain = Terrain(sources)

    az, d, vis, z0 = compute_coverage(terrain, a.lat, a.lon, a.ant_h, [f * FT for f in heights_ft], rmax,
                                      step, a.az_step, a.k, a.ground_m, a.agl)
    if getattr(terrain, "unknown", 0) > 0.02 and not a.flat:
        print(f"Attention : {terrain.unknown:.0%} des points sans donnée de relief (traités à 0 m).", file=sys.stderr)

    cont, vmax = write_kml(a.output, a.lat, a.lon, z0, a.ant_h, az, d, vis, heights_ft, a.k,
                           f"Couverture ADS-B/ADS-L {a.lat:.4f}, {a.lon:.4f}", a.agl, a.az_step, step)
    ref = "AGL" if a.agl else "AMSL"
    print(f"Sol récepteur : {z0:.0f} m ; antenne : {z0 + a.ant_h:.0f} m AMSL ; k = {a.k:.4f} ; niveaux {ref}")
    print(f"{'Niveau':>11} | portée continue km : {'min':>5} {'moy':>5} {'max':>5} | visible max | points visibles")
    for f, r, vm, v in zip(heights_ft, cont, vmax, vis):
        print(f"{f:>7g} ft  |{'':>21}{r.min() / 1000:5.1f} {r.mean() / 1000:5.1f} {r.max() / 1000:5.1f} | "
              f"{vm.max() / 1000:8.1f} km | {v.mean():8.1%}")
    print(f"KML : {a.output}")


if __name__ == "__main__":
    main()
