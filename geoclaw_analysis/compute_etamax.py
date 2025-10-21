"""Compute shoreline tsunami proxies (eta statistics) for all available runs."""

from __future__ import annotations

import argparse
import os
import re
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


DEFAULT_RESOLUTIONS = ["coarse", "fine"]
FGMAX_NUMBERS = (1, 2, 3)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute shoreline eta statistics from fgmax results."
    )
    parser.add_argument(
        "--resolutions",
        nargs="*",
        default=DEFAULT_RESOLUTIONS,
        help="Resolutions to process (default: coarse fine)",
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
        "--output-dir",
        default=None,
        help="Where to write summary files (default: alongside input data).",
    )
    return parser.parse_args()


def discover_run_sets(base_dir: str, resolutions: Sequence[str], magnitudes: Optional[Sequence[str]]) -> List[Tuple[str, str, str]]:
    combos: List[Tuple[str, str, str]] = []
    for entry in os.listdir(base_dir):
        match = re.match(r'^(?P<res>[A-Za-z0-9]+)_(?P<mw>[\d.]+)_runs$', entry)
        if not match:
            continue
        res = match.group('res')
        mw = match.group('mw')
        if resolutions and res not in resolutions:
            continue
        if magnitudes and mw not in magnitudes:
            continue
        combos.append((res, mw, entry))

    def sort_key(item: Tuple[str, str, str]):
        res, mw, _ = item
        try:
            return (res, float(mw))
        except ValueError:
            return (res, mw)

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
    base_dir = os.path.abspath(os.path.join(script_dir, '..', 'geoclaw_output', 'all_runs_npy_files'))
    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f"Expected directory not found: {base_dir}")

    explicit_runs = [int(r) for r in args.runs] if args.runs else None

    combos = discover_run_sets(base_dir, args.resolutions, args.magnitudes)
    if not combos:
        print("No matching run directories found.")
        return

    output_dir = args.output_dir or base_dir
    os.makedirs(output_dir, exist_ok=True)

    summary: Dict[Tuple[str, str, int], Dict[str, object]] = defaultdict(lambda: {
        'processed': [],
        'missing': []
    })

    for resolution, Mw, run_root_name in combos:
        run_root = os.path.join(base_dir, run_root_name)
        runs = discover_runs(run_root, explicit_runs)
        if not runs:
            print(f"[WARN] No runs discovered in {run_root}")
            continue

        for fgno in FGMAX_NUMBERS:
            b0_path = os.path.join(base_dir, f"{resolution}_xyB0_fg{fgno}.npy")
            if not os.path.exists(b0_path):
                print(f"[WARN] Missing B0 file {b0_path}; skipping {resolution}_{Mw} FG{fgno}.")
                continue

            xyB0 = np.load(b0_path)
            B0 = xyB0[:, 2]

            output_path = os.path.join(output_dir, f"{resolution}_{Mw}_etamax_fg{fgno}.txt")
            with open(output_path, 'w') as etamax_file:
                etamax_file.write('  run     min eta    mean eta    max eta (on shore)  hmax sum\n')

                for runno, run_dir_name in runs:
                    hmax_path = os.path.join(run_root, run_dir_name, f'hmax_fg{fgno}.npy')
                    if not os.path.exists(hmax_path):
                        print(f"[WARN] Missing hmax for {resolution}_{Mw} run {runno} FG{fgno}")
                        summary[(resolution, Mw, fgno)]['missing'].append(runno)
                        etamax_file.write(f"{runno:4d}    {np.nan:8.3f}    {np.nan:8.3f}    {np.nan:8.3f}    {np.nan:17.3f}\n")
                        continue

                    hmax = np.load(hmax_path)
                    if hmax.shape[0] != B0.shape[0]:
                        print(f"[WARN] Shape mismatch for {resolution}_{Mw} run {runno} FG{fgno}; skipping.")
                        summary[(resolution, Mw, fgno)]['missing'].append(runno)
                        etamax_file.write(f"{runno:4d}    {np.nan:8.3f}    {np.nan:8.3f}    {np.nan:8.3f}    {np.nan:17.3f}\n")
                        continue

                    emin, emean, emax, hmax_sum = compute_eta_metrics(B0, hmax)
                    summary[(resolution, Mw, fgno)]['processed'].append((runno, emin, emean, emax, hmax_sum))
                    etamax_file.write(
                        f"{runno:4d}    {emin:8.3f}    {emean:8.3f}    {emax:8.3f}    {hmax_sum:17.3f}\n"
                    )

            print(f"Created {output_path}")

    if not summary:
        print("No eta statistics were generated.")
        return

    print("\nSummary report:")
    for (resolution, Mw, fgno), data in sorted(summary.items()):
        processed = data['processed']
        missing = data['missing']
        label = f"{resolution}_{Mw} FG{fgno}"
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
