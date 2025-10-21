"""Generate coarse/fine scatter plots from fgmax shoreline metrics.

The script scans the ``all_runs_npy_files`` directory and only processes
cases where both coarse and fine data are present. Figures are written to the
``png_files`` directory (no interactive display). A summary of processed and
missing runs is printed at the end.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

HAS_MPL: Optional[bool] = None
plt = None
matplotlib = None


DEFAULT_RESOLUTIONS = ["coarse", "fine"]
FGMAX_NUMBERS = (1, 2, 3)
MW_COLORS = {8.6: 'm', 8.8: 'g', 9.0: 'b', 9.2: 'r'}

GLOBAL_OPTS: Dict[str, bool] = {}


def ensure_matplotlib() -> bool:
    global HAS_MPL, plt, matplotlib
    if HAS_MPL is not None:
        return HAS_MPL
    try:
        import matplotlib  # type: ignore
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt_module  # type: ignore
        plt = plt_module
        globals()['matplotlib'] = matplotlib
        HAS_MPL = True
    except Exception as err:  # pragma: no cover - optional dependency
        print(f"[WARN] Matplotlib unavailable: {err}")
        HAS_MPL = False
    return HAS_MPL


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Produce coarse/fine fgmax scatter plots from available data."
    )
    parser.add_argument(
        "--resolutions",
        nargs="*",
        default=DEFAULT_RESOLUTIONS,
        help="Resolutions to include (default: coarse fine).",
    )
    parser.add_argument(
        "--magnitudes",
        nargs="*",
        help="Optional list of magnitudes (e.g. 8.6 8.8) to include.",
    )
    parser.add_argument(
        "--runs",
        nargs="*",
        help="Optional explicit run numbers to include (e.g. 17 31).",
    )
    parser.add_argument(
        "--png-dir",
        default=None,
        help="Directory for saved figures (default: ../png_files).",
    )
    parser.add_argument(
        "--base-dir",
        default=None,
        help="Location of all_runs_npy_files (default: ../geoclaw_output/all_runs_npy_files).",
    )
    parser.add_argument(
        "--save-figures",
        action="store_true",
        help="Attempt to render PNG scatter plots in addition to data files.",
    )
    return parser.parse_args()


def discover_tables(base_dir: str, resolutions: Sequence[str], magnitudes: Optional[Sequence[str]]) -> Dict[Tuple[str, str, int], str]:
    pattern = os.path.join(base_dir, '*_etamax_fg?.txt')
    tables: Dict[Tuple[str, str, int], str] = {}
    for path in glob.iglob(pattern):
        entry = os.path.basename(path)
        match = re.match(r'^(?P<res>[A-Za-z0-9]+)_(?P<mw>[\d.]+)_etamax_fg(?P<fgno>[123])\.txt$', entry)
        if not match:
            continue
        res = match.group('res')
        mw = match.group('mw')
        fgno = int(match.group('fgno'))
        if resolutions and res not in resolutions:
            continue
        if magnitudes and mw not in magnitudes:
            continue
        tables[(res, mw, fgno)] = path
    return tables


def load_table(path: str) -> Dict[int, np.ndarray]:
    data = np.loadtxt(path, skiprows=1)
    if data.ndim == 1:
        data = data[None, :]
    table = {int(row[0]): row[1:] for row in data}
    return table


def ensure_output_dir(out_dir: str) -> str:
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


def save_scatter_data(out_dir: str, filename: str, data: Dict[float, Tuple[np.ndarray, np.ndarray]]) -> str:
    path = os.path.join(out_dir, filename)
    payload = {}
    for Mw, (xvals, yvals) in data.items():
        payload[f"Mw_{Mw:.1f}_x"] = xvals
        payload[f"Mw_{Mw:.1f}_y"] = yvals
    np.savez(path, **payload)
    return path + '.npz'


def save_scatter_plot(out_dir: str, filename: str, data: Dict[float, Tuple[np.ndarray, np.ndarray]], xlabel: str, ylabel: str, title: str) -> str:
    if not GLOBAL_OPTS.get('save_figures'):
        return ""
    if not ensure_matplotlib():
        return ""

    fig, ax = plt.subplots(figsize=(6, 6))
    qoi_min = -1.0
    qoi_max = -np.inf

    for Mw, (xvals, yvals) in sorted(data.items()):
        if xvals.size == 0:
            continue
        ax.scatter(xvals, yvals, color=MW_COLORS.get(Mw, 'k'), label=f'Mw = {Mw}')
        qoi_max = max(qoi_max, np.nanmax(xvals), np.nanmax(yvals))
        qoi_min = min(qoi_min, np.nanmin(xvals), np.nanmin(yvals))

    if qoi_max == -np.inf:
        plt.close(fig)
        return ""

    limits = [0.9 * qoi_min, 1.1 * qoi_max]
    ax.plot(limits, limits, 'k-')
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect('equal', 'box')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc='lower right')

    path = os.path.join(out_dir, f'{filename}.png')
    try:
        fig.savefig(path)
        return path
    except Exception as err:  # pragma: no cover - plotting fallback
        print(f"[WARN] Failed to save figure {path}: {err}")
        return ""
    finally:
        plt.close(fig)


def extract_qoi(table: Dict[int, np.ndarray], qoi: str) -> Dict[int, float]:
    column_map = {
        'etamin': 0,
        'etamean': 1,
        'etamax': 2,
        'hmax_sum': 3,
    }
    idx = column_map[qoi]
    return {run: row[idx] for run, row in table.items()}


def main() -> None:
    args = parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = args.base_dir or os.path.abspath(os.path.join(script_dir, '..', 'geoclaw_output', 'all_runs_npy_files'))
    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"Expected directory not found: {base_dir}")
    print(f"Scanning results in {base_dir}")

    GLOBAL_OPTS['save_figures'] = args.save_figures
    if GLOBAL_OPTS['save_figures']:
        ensure_matplotlib()

    explicit_runs = [int(r) for r in args.runs] if args.runs else None

    tables = discover_tables(base_dir, args.resolutions, args.magnitudes)
    if not tables:
        print("No etamax tables found.")
        return
    print(f"Found {len(tables)} etamax tables")

    output_dir = args.png_dir or os.path.abspath(os.path.join(script_dir, '..', 'png_files'))
    ensure_output_dir(output_dir)

    summary = defaultdict(lambda: {'processed': [], 'missing_coarse': [], 'missing_fine': []})

    coarse_fine_data = defaultdict(lambda: defaultdict(dict))  # (fgno,qoi)[Mw] = (coarse,fine)
    resolution_data = defaultdict(lambda: defaultdict(dict))    # (resolution,fgno)[Mw] = (x,y)

    tables_data = {}
    for key, path in tables.items():
        res, mw, fgno = key
        if res not in ('coarse', 'fine'):
            continue
        table = load_table(path)
        if explicit_runs is not None:
            table = {run: vals for run, vals in table.items() if run in explicit_runs}
        tables_data[key] = table

    mw_values = sorted({float(key[1]) for key in tables_data})

    for Mw in mw_values:
        Mw_str = f'{Mw}'
        for fgno in FGMAX_NUMBERS:
            coarse_table = tables_data.get(('coarse', Mw_str, fgno))
            fine_table = tables_data.get(('fine', Mw_str, fgno))
            if not coarse_table or not fine_table:
                continue

            common_runs = sorted(set(coarse_table) & set(fine_table))
            missing_coarse = sorted(set(fine_table) - set(coarse_table))
            missing_fine = sorted(set(coarse_table) - set(fine_table))

            summary[(Mw, fgno)]['missing_coarse'].extend(missing_coarse)
            summary[(Mw, fgno)]['missing_fine'].extend(missing_fine)

            if not common_runs:
                continue

            for qoi in ['etamin', 'etamean', 'etamax', 'hmax_sum']:
                coarse_vals = np.array([extract_qoi(coarse_table, qoi)[run] for run in common_runs])
                fine_vals = np.array([extract_qoi(fine_table, qoi)[run] for run in common_runs])
                coarse_fine_data[(fgno, qoi)][Mw] = (coarse_vals, fine_vals)

            for resolution in ['coarse', 'fine']:
                table = coarse_table if resolution == 'coarse' else fine_table
                etamean_vals = np.array([extract_qoi(table, 'etamean')[run] for run in common_runs])
                etamax_vals = np.array([extract_qoi(table, 'etamax')[run] for run in common_runs])
                resolution_data[(resolution, fgno)][Mw] = (etamean_vals, etamax_vals)
            summary[(Mw, fgno)]['processed'].append(('common_runs', '', len(common_runs)))

    data_files = []
    figure_files = []

    for (fgno, qoi), data in coarse_fine_data.items():
        data_path = save_scatter_data(output_dir, f'scatter_{qoi}_fg{fgno}', data)
        data_files.append(data_path)
        fig_path = save_scatter_plot(
            output_dir,
            f'scatter_{qoi}_fg{fgno}',
            data,
            'coarse',
            'fine',
            f'{qoi} over onshore points (fgmax {fgno})'
        )
        if fig_path:
            figure_files.append(fig_path)

    for (resolution, fgno), data in resolution_data.items():
        data_path = save_scatter_data(
            output_dir,
            f'scatter_etamean_etamax_{resolution}_fg{fgno}',
            data,
        )
        data_files.append(data_path)
        fig_path = save_scatter_plot(
            output_dir,
            f'scatter_etamean_etamax_{resolution}_fg{fgno}',
            data,
            f'etamean ({resolution})',
            f'etamax ({resolution})',
            f'etamax vs etamean ({resolution}) fgmax {fgno}'
        )
        if fig_path:
            figure_files.append(fig_path)

    if not summary:
        print("No matching coarse/fine pairs were processed.")
        return

    print("Summary report:")
    for (Mw, fgno), info in sorted(summary.items()):
        processed = info['processed']
        missing_coarse = sorted(set(info['missing_coarse']))
        missing_fine = sorted(set(info['missing_fine']))
        print(f"  Mw {Mw:.1f} fgmax {fgno}:")
        processed_runs = [entry[2] for entry in processed if entry[0] == 'common_runs']
        if processed_runs:
            print(f"    processed {processed_runs[0]} common run(s)")
        else:
            print("    no common runs processed")
        if missing_coarse:
            print(f"    missing coarse runs: {', '.join(str(r) for r in missing_coarse)}")
        if missing_fine:
            print(f"    missing fine runs: {', '.join(str(r) for r in missing_fine)}")
    if data_files:
        print("Saved scatter datasets:")
        for path in data_files:
            print(f"  {path}")
    if figure_files:
        print("Saved scatter figures:")
        for path in figure_files:
            print(f"  {path}")
    if not GLOBAL_OPTS.get('save_figures'):
        print("[INFO] Figures not requested; use --save-figures to generate PNGs.")
    elif not figure_files:
        print("[INFO] Matplotlib unavailable or no figures generated.")


if __name__ == '__main__':
    main()
