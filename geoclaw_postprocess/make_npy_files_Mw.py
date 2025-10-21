"""Generate per-run numpy summaries for GeoClaw fgmax (and optional gauges).

The script scans the existing ``geoclaw_output`` directory and processes any
``<resolution>_<Mw>/run_<id>`` subdirectories that contain fgmax output. It is
designed to work with both full production runs and small "test" batches where
only a handful of runs exist (e.g., there may be no ``run_0``).

Saved outputs are written under ``all_runs_npy_files/<resolution>_<Mw>_runs/``.
"""

from __future__ import annotations

import argparse
import glob
import importlib
import os
import re
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
from scipy.interpolate import interp1d

from clawpack.geoclaw import fgmax_tools


DEFAULT_RESOLUTIONS = ["coarse", "fine"]
FGMAX_NUMBERS = (1, 2, 3)
DTOTO_CANDIDATES = [
    ("dtopo.tt3", 3),
    ("dtopo1.tt3", 3),
    ("dtopo2.tt3", 3),
    ("dtopo.tt2", 2),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create numpy archives for GeoClaw fgmax (and optional gauges)."
    )
    parser.add_argument(
        "--resolutions",
        nargs="*",
        default=DEFAULT_RESOLUTIONS,
        help="Resolutions to process (default: coarse fine).",
    )
    parser.add_argument(
        "--magnitudes",
        nargs="*",
        help="Optional list of magnitudes (e.g. 8.6 8.8) to include.",
    )
    parser.add_argument(
        "--runs",
        nargs="*",
        help="Optional explicit run numbers to process (e.g. 17 31 76).",
    )
    parser.add_argument(
        "--make-gauge-npy",
        action="store_true",
        help="Also create gauge numpy files (disabled by default).",
    )
    return parser.parse_args()


def parse_fgmax_legacy(run_fg: fgmax_tools.FGmaxGrid, fgno: int, resolution: str) -> None:
    legacy_path = os.path.join('..', 'geoclaw_driver', f'fgmax{fgno}_{resolution}.txt')
    if not os.path.exists(legacy_path):
        raise FileNotFoundError(f"FGmax definition missing: {legacy_path}")

    tokens: List[str] = []
    detected_fgno: Optional[int] = None
    with open(legacy_path) as handle:
        for raw in handle:
            main, _, comment = raw.partition('#')
            main = main.strip()
            comment = comment.strip().lower()
            if 'fgno' in comment and detected_fgno is None:
                if main:
                    detected_fgno = int(float(main))
                continue
            if main:
                tokens.append(main)

    if len(tokens) < 7:
        raise ValueError(f"Unexpected format in {legacy_path}")

    idx = 0
    run_fg.fgno = detected_fgno if detected_fgno is not None else fgno
    run_fg.tstart_max = float(tokens[idx]); idx += 1
    run_fg.tend_max = float(tokens[idx]); idx += 1
    run_fg.dt_check = float(tokens[idx]); idx += 1
    run_fg.min_level_check = int(float(tokens[idx])); idx += 1
    run_fg.arrival_tol = float(tokens[idx]); idx += 1
    run_fg.interp_method = int(float(tokens[idx])); idx += 1
    run_fg.point_style = point_style = int(float(tokens[idx])); idx += 1

    if point_style == 0:
        run_fg.npts = int(float(tokens[idx])); idx += 1
        if run_fg.npts == 0:
            run_fg.xy_fname = tokens[idx].strip('"')
        else:
            coords = [list(map(float, tok.split())) for tok in tokens[idx:idx + run_fg.npts]]
            run_fg.X = np.array([xy[0] for xy in coords])
            run_fg.Y = np.array([xy[1] for xy in coords])
    elif point_style == 1:
        run_fg.npts = int(float(tokens[idx])); idx += 1
        run_fg.x1, run_fg.y1 = map(float, tokens[idx].split()); idx += 1
        run_fg.x2, run_fg.y2 = map(float, tokens[idx].split()); idx += 1
    elif point_style == 2:
        run_fg.nx, run_fg.ny = [int(float(v)) for v in tokens[idx].split()]; idx += 1
        run_fg.npts = run_fg.nx * run_fg.ny
        run_fg.x1, run_fg.y1 = map(float, tokens[idx].split()); idx += 1
        run_fg.x2, run_fg.y2 = map(float, tokens[idx].split()); idx += 1
    elif point_style == 3:
        run_fg.n12, run_fg.n23 = [int(float(v)) for v in tokens[idx].split()]; idx += 1
        run_fg.npts = run_fg.n12 * run_fg.n23
        run_fg.x1, run_fg.y1 = map(float, tokens[idx].split()); idx += 1
        run_fg.x2, run_fg.y2 = map(float, tokens[idx].split()); idx += 1
        run_fg.x3, run_fg.y3 = map(float, tokens[idx].split()); idx += 1
        run_fg.x4, run_fg.y4 = map(float, tokens[idx].split()); idx += 1
    else:
        raise ValueError(f"Point style {point_style} not supported in {legacy_path}")


def load_fgmax_grid(run_fg: fgmax_tools.FGmaxGrid, fgno: int, resolution: str) -> None:
    grids_file = os.path.join('..', 'geoclaw_driver', f'fgmax_grids_{resolution}.data')
    if os.path.exists(grids_file):
        try:
            run_fg.read_fgmax_grids_data(fgno, data_file=grids_file)
            return
        except ValueError:
            pass
    parse_fgmax_legacy(run_fg, fgno, resolution)


def load_dtopo(rundir: str):
    from clawpack.geoclaw.dtopotools import DTopography

    for fname, topo_type in DTOTO_CANDIDATES:
        path = os.path.join(rundir, fname)
        if os.path.exists(path):
            dtopo = DTopography()
            dtopo.read(path, dtopo_type=topo_type)
            return dtopo
    raise FileNotFoundError(f"No dtopo file found in {rundir}")


def discover_magnitude_dirs(resolution: str, magnitude_filter: Optional[Sequence[str]]) -> List[Tuple[str, str]]:
    entries = []
    for entry in os.listdir('.'):
        if not entry.startswith(f"{resolution}_"):
            continue
        suffix = entry[len(resolution) + 1:]
        if suffix.upper() == 'B0':
            continue
        if magnitude_filter and suffix not in magnitude_filter:
            continue
        if os.path.isdir(entry):
            entries.append((suffix, entry))

    def mag_key(item: Tuple[str, str]):
        try:
            return float(item[0])
        except ValueError:
            return item[0]

    entries.sort(key=mag_key)
    return entries


def discover_run_dirs(mw_dir: str, explicit_runs: Optional[Sequence[int]]) -> List[str]:
    run_dirs = []
    for entry in os.listdir(mw_dir):
        full = os.path.join(mw_dir, entry)
        if os.path.isdir(full) and re.match(r'^run_\d+$', entry):
            run_dirs.append(entry)

    if explicit_runs is not None:
        allowed = {int(r) for r in explicit_runs}
        run_dirs = [d for d in run_dirs if int(d.split('_')[1]) in allowed]

    run_dirs.sort(key=lambda name: int(name.split('_')[1]))
    return run_dirs


def main() -> None:
    args = parse_args()

    magnitude_filter = set(args.magnitudes) if args.magnitudes else None
    explicit_runs = [int(r) for r in args.runs] if args.runs else None
    gauges_module = None

    script_dir = os.path.dirname(os.path.abspath(__file__))
    geoclaw_output = os.path.abspath(os.path.join(script_dir, '..', 'geoclaw_output'))

    if not os.path.isdir(geoclaw_output):
        raise FileNotFoundError(f"Expected directory not found: {geoclaw_output}")

    os.chdir(geoclaw_output)
    print(f"Working in directory: {geoclaw_output}")

    all_runs_dir = 'all_runs_npy_files'
    os.makedirs(all_runs_dir, exist_ok=True)

    # Crescent City location:
    xcc = -124.1838
    ycc = 41.7456

    for resolution in args.resolutions:
        magnitude_dirs = discover_magnitude_dirs(resolution, magnitude_filter)
        if not magnitude_dirs:
            print(f"[WARN] No magnitude directories found for resolution '{resolution}'.")
            continue

        for Mw, Mw_dir in magnitude_dirs:
            run_dirs = discover_run_dirs(Mw_dir, explicit_runs)
            if not run_dirs:
                print(f"[WARN] No run directories found in {Mw_dir}.")
                continue

            dtopo_source = None
            for rundir in run_dirs:
                candidate = os.path.join(Mw_dir, rundir)
                try:
                    load_dtopo(candidate)
                except FileNotFoundError:
                    continue
                else:
                    dtopo_source = candidate
                    break

            if dtopo_source is None:
                print(f"[WARN] No dtopo files found for {Mw_dir}; skipping.")
                continue

            print(f"Using dtopo files from {dtopo_source}")
            dtopo_ref = load_dtopo(dtopo_source)
            i1cc = np.where(dtopo_ref.x < xcc)[0].max()
            j1cc = np.where(dtopo_ref.y < ycc)[0].max()
            a1cc = (xcc - dtopo_ref.x[i1cc]) / (dtopo_ref.x[i1cc + 1] - dtopo_ref.x[i1cc])
            a2cc = (ycc - dtopo_ref.y[j1cc]) / (dtopo_ref.y[j1cc + 1] - dtopo_ref.y[j1cc])
            if (a1cc < 0.) or (a1cc > 1.) or (a2cc < 0.) or (a2cc > 1.):
                print('*** Interpolation to CC not correct!')

            def dZ_CrescentCity(dtopo):
                dZr = dtopo.dZ[0, :, :]
                dzy1 = (1. - a1cc) * dZr[j1cc, i1cc] + a1cc * dZr[j1cc, i1cc + 1]
                dzy2 = (1. - a1cc) * dZr[j1cc + 1, i1cc] + a1cc * dZr[j1cc + 1, i1cc + 1]
                dzcc = (1. - a2cc) * dzy2 + a2cc * dzy1
                return dzcc

            new_runs_dir = os.path.join(all_runs_dir, f'{resolution}_{Mw}_runs')
            os.makedirs(new_runs_dir, exist_ok=True)

            with open(os.path.join(new_runs_dir, 'dzCC.txt'), 'w') as dzfile:
                for rundir in run_dirs:
                    runno = int(rundir.split('_')[1])
                    fullrundir = os.path.join(Mw_dir, rundir)
                    outdir = os.path.join(fullrundir, '_output')
                    if not os.path.isdir(outdir):
                        print(f"[WARN] Missing _output directory for {fullrundir}; skipping.")
                        continue
                    if not glob.glob(os.path.join(outdir, 'fgmax000*.txt')):
                        print(f"[WARN] Missing fgmax outputs in {fullrundir}; skipping.")
                        continue

                    newrundir = os.path.join(new_runs_dir, rundir)
                    os.makedirs(newrundir, exist_ok=True)
                    print(f"Processing {fullrundir}")

                    try:
                        dtopo = load_dtopo(fullrundir)
                    except FileNotFoundError as err:
                        print(f"[WARN] {err}; skipping run {rundir}.")
                        continue

                    dzCC = dZ_CrescentCity(dtopo)
                    dzfile.write(f"{runno:4d}   {dzCC:13.8f}\n")
                    dz_path = os.path.join(newrundir, 'dzCC.txt')
                    np.savetxt(dz_path, np.array([dzCC]))
                    print("Created", dz_path)

                    for fgno in FGMAX_NUMBERS:
                        fg = fgmax_tools.FGmaxGrid()
                        load_fgmax_grid(fg, fgno, resolution)
                        fg.read_output(fgno, outdir)
                        hmax = np.reshape(fg.h, -1, order='F')
                        fg_path = os.path.join(newrundir, f'hmax_fg{fgno}.npy')
                        np.save(fg_path, np.array(hmax))
                        print("Created", fg_path)

                    if args.make_gauge_npy:
                        if gauges_module is None:
                            gauges_module = importlib.import_module('clawpack.pyclaw.gauges')
                        gauge_files = glob.glob(os.path.join(outdir, 'gauge0*txt'))
                        gaugenos = sorted(int(os.path.basename(f)[-9:-4]) for f in gauge_files)
                        if runno == 0 and gaugenos:
                            print("Found gauges", gaugenos)

                        tt = np.linspace(0, 8900, 891)
                        for gaugeno in gaugenos:
                            g = gauges_module.GaugeSolution(gaugeno, outdir)
                            eta_func = interp1d(g.t, g.q[3, :], bounds_error=False)
                            eta_t = eta_func(tt)
                            h_func = interp1d(g.t, g.q[0, :], bounds_error=False)
                            h_t = h_func(tt)
                            h_eta = np.vstack((h_t, eta_t)).T
                            gauge_path = os.path.join(newrundir, f'gauge{str(gaugeno).zfill(5)}.npy')
                            np.save(gauge_path, h_eta)
                            print(f"Created {gauge_path}")
                        if gaugenos:
                            print("Created gauge npy files")


if __name__ == '__main__':
    main()
