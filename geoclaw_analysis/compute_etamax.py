"""Compute shoreline tsunami proxies (eta statistics) for all available runs."""

from __future__ import annotations

import argparse
import functools
import os
import re
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from clawpack.geoclaw import fgmax_tools

FGMAX_NUMBERS = (1, 2, 3)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute shoreline eta statistics from fgmax results."
    )
    parser.add_argument(
        "--magnitudes",
        nargs="*",
        help="Optional list of moment magnitudes (e.g. 8.0 8.6 9.0) to include.",
    )
    parser.add_argument(
        "--runs",
        nargs="*",
        help="Optional explicit run numbers to include (e.g. 17 31).",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Where to write summary files (default: alongside input data).",
    )
    return parser.parse_args()


def discover_run_sets(base_dir: str, magnitudes: Optional[Sequence[str]]) -> List[Tuple[str, str, str]]:
    combos: List[Tuple[str, str, str]] = []
    for entry in os.listdir(base_dir):
        if not entry.startswith("Mw_") or not entry.endswith("_runs"):
            continue
        suffix = entry[3:-5]  # strip 'Mw_' and '_runs'
        mw_label = _decode_mw_suffix(suffix)
        if magnitudes and mw_label not in magnitudes:
            continue
        raw_dir = f"Mw_{suffix}"
        combos.append((mw_label, entry, raw_dir))

    def sort_key(item: Tuple[str, str, str]):
        mw_label, _, _ = item
        try:
            return float(mw_label)
        except ValueError:
            return mw_label

    combos.sort(key=sort_key)
    return combos


def discover_runs(run_root: str, explicit_runs: Optional[Sequence[int]]) -> List[Tuple[int, str]]:
    run_entries: List[Tuple[int, str]] = []
    for entry in os.listdir(run_root):
        match = re.match(r'^run_(\d+)$', entry)
        if match:
            run_no = int(match.group(1))
            if explicit_runs is not None and run_no not in explicit_runs:
                continue
            run_entries.append((run_no, entry))
    run_entries.sort(key=lambda item: item[0])
    return run_entries


def _decode_mw_suffix(suffix: str) -> str:
    try:
        value = float(suffix) / 10.0
        return f"{value:.1f}"
    except ValueError:
        return suffix


@functools.lru_cache(maxsize=None)
def load_b0(raw_root: str, fgno: int) -> Optional[np.ndarray]:
    """Load and flatten the fgmax topography (B0) array for the given magnitude and grid."""
    if not os.path.isdir(raw_root):
        return None

    run_candidates: List[str] = [
        entry
        for entry in os.listdir(raw_root)
        if re.match(r'^run_\d+$', entry)
    ]
    run_candidates.sort(key=lambda name: int(name.split('_')[1]))

    for run_name in run_candidates:
        grids_file = os.path.join(raw_root, run_name, "fgmax_grids.data")
        out_dir = os.path.join(raw_root, run_name, "_output")
        fg_file = os.path.join(out_dir, f"fgmax{fgno:04d}.txt")
        if not (os.path.exists(grids_file) and os.path.isdir(out_dir) and os.path.exists(fg_file)):
            continue

        fg = fgmax_tools.FGmaxGrid()
        try:
            fg.read_fgmax_grids_data(fgno, data_file=grids_file)
            fg.read_output(outdir=out_dir)
        except (IOError, ValueError):
            continue

        B = np.ma.filled(fg.B, np.nan)
        return B.reshape(-1, order='F')

    return None


def compute_eta_metrics(B0: np.ndarray, hmax: np.ndarray) -> Tuple[float, float, float, float]:
    eta = np.where(hmax > 0.0, B0 + hmax, np.nan)
    eta_shore = np.where((hmax > 0.0) & (B0 > 0.0), eta, np.nan)

    valid = np.isfinite(eta_shore)
    if not np.any(valid):
        emin = emean = emax = 0.0
    else:
        emin = float(np.nanmin(eta_shore))
        emax = float(np.nanmax(eta_shore))
        emean = float(np.nanmean(eta_shore))

    hmax_sum = float(np.sum(np.where(B0 > 0.0, hmax, 0.0)))
    return emin, emean, emax, hmax_sum


def main() -> None:
    args = parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_root = os.path.abspath(os.path.join(script_dir, '..', 'geoclaw_output'))
    base_dir = os.path.join(output_root, 'all_runs_npy_files')
    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"Expected directory not found: {base_dir}")

    explicit_runs = [int(r) for r in args.runs] if args.runs else None

    combos = discover_run_sets(base_dir, args.magnitudes)
    if not combos:
        print("No matching run directories found.")
        return

    output_dir = args.output_dir or base_dir
    os.makedirs(output_dir, exist_ok=True)

    summary: Dict[Tuple[str, int], Dict[str, object]] = defaultdict(lambda: {
        'processed': [],
        'missing': []
    })

    for Mw, run_root_name, raw_dir_name in combos:
        run_root = os.path.join(base_dir, run_root_name)
        runs = discover_runs(run_root, explicit_runs)
        if not runs:
            print(f"[WARN] No runs discovered in {run_root}")
            continue

        raw_fg_root = os.path.join(output_root, raw_dir_name)

        for fgno in FGMAX_NUMBERS:
            B0 = load_b0(raw_fg_root, fgno)
            if B0 is None:
                print(f"[WARN] Unable to load B0 for {raw_dir_name} FG{fgno}; skipping.")
                continue

            output_path = os.path.join(output_dir, f"{run_root_name}_etamax_fg{fgno}.txt")
            with open(output_path, 'w') as etamax_file:
                etamax_file.write('  run     min eta    mean eta    max eta (on shore)  hmax sum\n')

                for runno, run_dir_name in runs:
                    hmax_path = os.path.join(run_root, run_dir_name, f'hmax_fg{fgno}.npy')
                    if not os.path.exists(hmax_path):
                        print(f"[WARN] Missing hmax for {raw_dir_name} run {runno} FG{fgno}")
                        summary[(Mw, fgno)]['missing'].append(runno)
                        etamax_file.write(f"{runno:4d}    {np.nan:8.3f}    {np.nan:8.3f}    {np.nan:8.3f}    {np.nan:17.3f}\n")
                        continue

                    hmax = np.load(hmax_path)
                    if hmax.shape[0] != B0.shape[0]:
                        print(f"[WARN] Shape mismatch for {raw_dir_name} run {runno} FG{fgno}; skipping.")
                        summary[(Mw, fgno)]['missing'].append(runno)
                        etamax_file.write(f"{runno:4d}    {np.nan:8.3f}    {np.nan:8.3f}    {np.nan:8.3f}    {np.nan:17.3f}\n")
                        continue

                    emin, emean, emax, hmax_sum = compute_eta_metrics(B0, hmax)
                    summary[(Mw, fgno)]['processed'].append((runno, emin, emean, emax, hmax_sum))
                    etamax_file.write(
                        f"{runno:4d}    {emin:8.3f}    {emean:8.3f}    {emax:8.3f}    {hmax_sum:17.3f}\n"
                    )

            print(f"Created {output_path}")

    if not summary:
        print("No eta statistics were generated.")
        return

    print("\nSummary report:")
    for (Mw, fgno), data in sorted(summary.items()):
        processed = data['processed']
        missing = data['missing']
        label = f"Mw {Mw} FG{fgno}"
        print(f"  {label}: processed {len(processed)} run(s)")
        if processed:
            values = np.array([p[1:] for p in processed])
            means = values.mean(axis=0)
            mins = values.min(axis=0)
            maxs = values.max(axis=0)
            fields = ['min eta', 'mean eta', 'max eta', 'hmax sum']
            for field, mean_v, min_v, max_v in zip(fields, means, mins, maxs):
                print(f"    {field:<10s} mean={mean_v:9.3f} min={min_v:9.3f} max={max_v:9.3f}")
        if missing:
            print(f"    Missing runs: {', '.join(str(r) for r in missing)}")


if __name__ == '__main__':
    main()
