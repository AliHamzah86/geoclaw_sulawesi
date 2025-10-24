"""Compute tsunami source proxies for all available GeoClaw runs.

For every ``Mw_<code>/run_<id>`` directory discovered beneath
``../geoclaw_output`` the script evaluates:

* seafloor displacement at Crescent City (dzCC)
* maximum offshore uplift/subsidence (eta_max)
* displaced water volume (km^3)
* potential energy (peta-joules)

It supports Clawpack 5.13 layouts, automatically locates dtopo files,
and tolerates partial test batches where only a handful of runs exist.

Results are written to ``geoclaw_output/all_runs_npy_files/dtopo_proxies_Mw_<code>.txt``
and a concise summary is printed at the end.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import numpy.ma as ma
from scipy.interpolate import RectBivariateSpline

from clawpack.geoclaw import topotools


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
RUN_CONFIG_PATH = os.path.join(PROJECT_ROOT, 'run_config.json')
DEFAULT_DATA_DIR = os.path.join(PROJECT_ROOT, 'DataFiles')
DEFAULT_GEOCLA_OUTPUT = os.path.join(PROJECT_ROOT, 'geoclaw_output')


DTOTO_CANDIDATES = [
    ("dtopo.tt3", 3),
    ("dtopo1.tt3", 3),
    ("dtopo2.tt3", 3),
    ("dtopo.tt2", 2),
]
SUMMARY_HEADER = (
    "run", "dzCC_m", "eta_max_m", "PE_pJ", "dVolume_km3"
)


def _resolve_defaults() -> Tuple[str, str, str, Dict[str, object]]:
    data_dir = DEFAULT_DATA_DIR
    output_base = DEFAULT_GEOCLA_OUTPUT
    config: Dict[str, object] = {}

    try:
        with open(RUN_CONFIG_PATH, 'r', encoding='utf-8') as handle:
            config = json.load(handle)
            if config.get('data_dir'):
                data_dir = os.path.normpath(os.path.expanduser(config['data_dir']))
            if config.get('output_base'):
                output_base = os.path.normpath(os.path.expanduser(config['output_base']))
    except FileNotFoundError:
        config = {}

    topo_path: Optional[str] = None
    topography_entries = config.get('topography') if isinstance(config.get('topography'), list) else []
    fallback_name = 'etopo1_-130_-124_38_45_1min.asc'

    if isinstance(topography_entries, list):
        for entry in topography_entries:
            if not isinstance(entry, dict):
                continue
            candidate = entry.get('path')
            if not candidate:
                continue
            candidate_expanded = os.path.expanduser(candidate)
            if not os.path.isabs(candidate_expanded):
                candidate_expanded = os.path.normpath(os.path.join(data_dir, candidate_expanded))
            if os.path.exists(candidate_expanded):
                topo_path = candidate_expanded
                break

            alt = os.path.normpath(os.path.join(DEFAULT_DATA_DIR, os.path.basename(candidate_expanded)))
            if os.path.exists(alt):
                topo_path = alt
                break

    if not topo_path:
        fallback = os.path.join(data_dir, fallback_name)
        if not os.path.exists(fallback):
            fallback = os.path.join(DEFAULT_DATA_DIR, fallback_name)
        topo_path = os.path.normpath(fallback)

    return topo_path, data_dir, output_base, config


DEFAULT_TOPO_PATH, CONFIG_DATA_DIR, CONFIG_OUTPUT_BASE, _CONFIG_CONTENT = _resolve_defaults()
DEFAULT_PROXIES_DIR = os.path.join(CONFIG_OUTPUT_BASE, 'all_runs_npy_files')


def resolve_topography_path(path: str) -> str:
    expanded = os.path.expanduser(path)
    if os.path.isabs(expanded):
        return os.path.normpath(expanded)

    search_roots = [
        CONFIG_DATA_DIR,
        DEFAULT_DATA_DIR,
        PROJECT_ROOT,
    ]

    for root in search_roots:
        if not root:
            continue
        candidate = os.path.normpath(os.path.join(root, expanded))
        if os.path.exists(candidate):
            return candidate

    fallback_root = CONFIG_DATA_DIR or DEFAULT_DATA_DIR
    return os.path.normpath(os.path.join(fallback_root, expanded))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute dtopo-based tsunami source proxies for GeoClaw runs."
    )
    parser.add_argument(
        "--magnitudes",
        nargs="*",
        help="Optional list of magnitudes (e.g. 8.6 8.8) to include.",
    )
    parser.add_argument(
        "--runs",
        nargs="*",
        help="Optional explicit run numbers to include (e.g. 17 31 76).",
    )
    parser.add_argument(
        "--topo-file",
        default=DEFAULT_TOPO_PATH,
        help="Topography file for masking (default: %(default)s).",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_PROXIES_DIR,
        help="Output directory for proxy tables (default: %(default)s).",
    )
    parser.add_argument(
        "--crescent-city",
        nargs=2,
        type=float,
        default=(-124.1838, 41.7456),
        metavar=("LON", "LAT"),
        help="Longitude/latitude of Crescent City reference point.",
    )
    parser.add_argument(
        "--info",
        action="store_true",
        help="Show resolved defaults (topography, output directories) and exit.",
    )
    return parser.parse_args()


def discover_magnitude_dirs(magnitude_filter: Optional[Sequence[str]]) -> List[Tuple[str, str]]:
    candidates: List[Tuple[str, str]] = []
    for entry in os.listdir('.'):
        if not entry.startswith("Mw_"):
            continue
        suffix = entry[3:]
        Mw_value = _decode_mw_suffix(suffix)
        if magnitude_filter and Mw_value not in magnitude_filter:
            continue
        if os.path.isdir(entry):
            candidates.append((Mw_value, entry))

    def magnitude_key(item: Tuple[str, str]):
        try:
            return float(item[0])
        except ValueError:
            return item[0]

    candidates.sort(key=magnitude_key)
    return candidates


def discover_run_dirs(mw_dir: str, explicit_runs: Optional[Sequence[int]]) -> List[str]:
    runs = []
    for entry in os.listdir(mw_dir):
        full = os.path.join(mw_dir, entry)
        if os.path.isdir(full) and re.match(r'^run_\d+$', entry):
            runs.append(entry)

    if explicit_runs is not None:
        allowed = {int(r) for r in explicit_runs}
        runs = [r for r in runs if int(r.split('_')[1]) in allowed]

    runs.sort(key=lambda name: int(name.split('_')[1]))
    return runs


def _decode_mw_suffix(suffix: str) -> str:
    """Convert directory suffix (e.g. '80') into magnitude string '8.0'."""
    try:
        value = float(suffix) / 10.0
        return f"{value:.1f}"
    except ValueError:
        return suffix


def find_dtopo(rundir: str):
    from clawpack.geoclaw.dtopotools import DTopography

    for fname, topo_type in DTOTO_CANDIDATES:
        path = os.path.join(rundir, fname)
        if os.path.exists(path):
            dtopo = DTopography()
            dtopo.read(path, dtopo_type=topo_type)
            return dtopo
    raise FileNotFoundError(f"No dtopo file found in {rundir}")


def prepare_topography(topo_path: str) -> Tuple[RectBivariateSpline, np.ndarray]:
    topo = topotools.Topography()
    topo.read(topo_path, 3)
    spline = RectBivariateSpline(topo.y, topo.x, topo.Z)
    return spline, topo.Z


def compute_cell_areas(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, float]:
    if len(x) < 2 or len(y) < 2:
        raise ValueError("Insufficient dtopo grid resolution to compute areas")

    deg2m = 111e3
    dy = abs(y[1] - y[0]) * deg2m
    X, Y = np.meshgrid(x, y)
    dx = abs(x[1] - x[0]) * deg2m * np.cos(np.deg2rad(Y))
    return dx * dy, dy


def interpolate_dz_at_location(dtopo, lon: float, lat: float) -> Optional[float]:
    x = dtopo.x
    y = dtopo.y
    if not (x[0] <= lon <= x[-1]) or not (y[0] <= lat <= y[-1]):
        return None

    i1 = np.searchsorted(x, lon) - 1
    j1 = np.searchsorted(y, lat) - 1
    i1 = np.clip(i1, 0, len(x) - 2)
    j1 = np.clip(j1, 0, len(y) - 2)

    a1 = (lon - x[i1]) / (x[i1 + 1] - x[i1])
    a2 = (lat - y[j1]) / (y[j1 + 1] - y[j1])

    dZr = dtopo.dZ[0, :, :]
    dzy1 = (1. - a1) * dZr[j1, i1] + a1 * dZr[j1, i1 + 1]
    dzy2 = (1. - a1) * dZr[j1 + 1, i1] + a1 * dZr[j1 + 1, i1 + 1]
    return (1. - a2) * dzy2 + a2 * dzy1


def compute_proxies(
    dtopo,
    topo_spline: RectBivariateSpline,
    cell_area: np.ndarray,
    lon: float,
    lat: float,
) -> Optional[Tuple[float, float, float, float]]:
    dzcc = interpolate_dz_at_location(dtopo, lon, lat)
    if dzcc is None:
        return None

    dZr = dtopo.dZ[0, :, :]

    topo_on_grid = topo_spline(dtopo.y, dtopo.x)
    water_mask = topo_on_grid > 0  # mask land
    eta = ma.masked_where(water_mask, dZr)

    area = ma.masked_where(water_mask, cell_area)

    grav = 9.81
    rho_water = 1000.0

    energy = (eta ** 2 * area).sum() * grav * rho_water * 1e-15
    volume = (eta * area).sum() * 1e-9
    eta_max = eta.max()

    return float(dzcc), float(eta_max), float(energy), float(volume)


def main() -> None:
    args = parse_args()

    magnitude_filter = set(args.magnitudes) if args.magnitudes else None
    explicit_runs = [int(r) for r in args.runs] if args.runs else None

    geoclaw_output = CONFIG_OUTPUT_BASE
    topo_path = resolve_topography_path(args.topo_file)
    output_base = os.path.expanduser(args.output_dir)
    if not os.path.isabs(output_base):
        output_base = os.path.normpath(os.path.join(CONFIG_OUTPUT_BASE, output_base))

    if args.info:
        print("Configuration summary")
        print(f"  run_config.json : {RUN_CONFIG_PATH} ({'exists' if os.path.exists(RUN_CONFIG_PATH) else 'missing'})")
        print(f"  data directory  : {CONFIG_DATA_DIR}")
        print(f"  GeoClaw output  : {geoclaw_output}")
        print(f"  proxies output  : {output_base}")
        print(f"  topography file : {topo_path} ({'exists' if os.path.exists(topo_path) else 'missing'})")
        lon_info, lat_info = args.crescent_city
        print(f"  Crescent City   : lon={lon_info} lat={lat_info}")
        return

    if not os.path.isdir(geoclaw_output):
        raise FileNotFoundError(f"Expected directory not found: {geoclaw_output}")
    if not os.path.exists(topo_path):
        raise FileNotFoundError(f"Topography file not found: {topo_path}")

    topo_spline, _ = prepare_topography(topo_path)

    os.chdir(geoclaw_output)
    print(f"Working in directory: {geoclaw_output}")

    os.makedirs(output_base, exist_ok=True)

    lon_cc, lat_cc = args.crescent_city

    summary: Dict[str, List[Tuple[int, float, float, float, float]]] = defaultdict(list)

    magnitude_dirs = discover_magnitude_dirs(magnitude_filter)
    if not magnitude_dirs:
        print("[WARN] No magnitude directories found.")
        return

    for Mw, Mw_dir in magnitude_dirs:
        run_dirs = discover_run_dirs(Mw_dir, explicit_runs)
        if not run_dirs:
            print(f"[WARN] No run directories found in {Mw_dir}.")
            continue

        out_path = os.path.join(output_base, f'dtopo_proxies_{Mw_dir}.txt')
        if os.path.exists(out_path):
            os.remove(out_path)
            print(f"Removed existing file: {out_path}")
        with open(out_path, 'w') as proxyfile:
            proxyfile.write(
                "  run       dzCC (m)      eta_max (m)      PE (pJ)      dVolume (km^3)\n"
            )

            for rundir in run_dirs:
                runno = int(rundir.split('_')[1])
                full_rundir = os.path.join(Mw_dir, rundir)
                try:
                    dtopo = find_dtopo(full_rundir)
                except FileNotFoundError:
                    print(f"[WARN] No dtopo file for {full_rundir}; skipping.")
                    continue

                cell_area, _ = compute_cell_areas(dtopo.x, dtopo.y)

                proxies = compute_proxies(dtopo, topo_spline, cell_area, lon_cc, lat_cc)
                if proxies is None:
                    print(
                        f"[WARN] Crescent City ({lon_cc}, {lat_cc}) outside domain for {full_rundir}; skipping."
                    )
                    continue

                dzcc, eta_max, energy, volume = proxies
                proxyfile.write(
                    f"{runno:4d}   {dzcc:13.8f}  {eta_max:13.8f}  {energy:13.8f}  {volume:13.8f}\n"
                )
                summary[Mw].append((runno, dzcc, eta_max, energy, volume))

        print(f"Created {out_path}")

    if not summary:
        print("No dtopo proxies were generated.")
        return

    print("\nSummary report:")
    for Mw, rows in sorted(summary.items()):
        values = np.array([row[1:] for row in rows])  # exclude run number
        counts = values.shape[0]
        means = values.mean(axis=0)
        mins = values.min(axis=0)
        maxs = values.max(axis=0)
        print(f"  Mw {Mw}: {counts} run(s)")
        for label, mean, min_v, max_v in zip(SUMMARY_HEADER[1:], means, mins, maxs):
            print(f"    {label:11s} mean={mean:9.4f} min={min_v:9.4f} max={max_v:9.4f}")


if __name__ == '__main__':
    main()
