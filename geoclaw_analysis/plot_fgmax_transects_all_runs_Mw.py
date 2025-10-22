"""
Plot transect cross sections for all available GeoClaw runs.

This script traverses the ``geoclaw_output`` directory, gathering fgmax
results for configured Mw scenarios and plotting shoreline transects for
fgmax grids 1 and 2 across a subset of runs.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from clawpack.geoclaw import fgmax_tools  # noqa: E402

BASE_DIR = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = BASE_DIR / "geoclaw_output"
ALL_RUNS_DIR = OUTPUT_ROOT / "all_runs_npy_files"
PNG_ROOT = BASE_DIR / "png_files"
PNG_ROOT.mkdir(parents=True, exist_ok=True)

TRANSECT_FGMAX = (1, 2)
MAX_RUNS_PER_PLOT = 50
ETA_COLOR = "r"


def decode_mw(dirname: str) -> Optional[Tuple[str, Path]]:
    if not dirname.startswith("Mw_") or not dirname.endswith("_runs"):
        return None
    suffix = dirname[3:-5]
    try:
        value = float(suffix) / 10.0
        label = f"{value:.1f}"
    except ValueError:
        label = suffix
    return label, Path(dirname)


def discover_mw_runs() -> List[Tuple[str, Path]]:
    entries: List[Tuple[str, Path]] = []
    if not ALL_RUNS_DIR.exists():
        return entries
    for child in sorted(ALL_RUNS_DIR.iterdir()):
        decoded = decode_mw(child.name)
        if decoded is not None and child.is_dir():
            entries.append((decoded[0], child))
    entries.sort(key=lambda item: float(item[0]) if item[0].replace(".", "", 1).isdigit() else item[0])
    return entries


def iter_run_dirs(run_root: Path) -> Iterable[Tuple[int, Path]]:
    for child in sorted(run_root.iterdir()):
        if child.is_dir() and child.name.startswith("run_"):
            try:
                run_number = int(child.name.split("_")[1])
            except (IndexError, ValueError):
                continue
            yield run_number, child


def load_fgmax_geometry(raw_root: Path, fgno: int) -> Optional[Dict[str, np.ndarray]]:
    if not raw_root.is_dir():
        return None
    for _, run_dir in iter_run_dirs(raw_root):
        grids_file = run_dir / "fgmax_grids.data"
        out_dir = run_dir / "_output"
        fg_file = out_dir / f"fgmax{fgno:04d}.txt"
        if not (grids_file.exists() and out_dir.is_dir() and fg_file.exists()):
            continue
        fg = fgmax_tools.FGmaxGrid()
        try:
            fg.read_fgmax_grids_data(fgno, data_file=str(grids_file))
            fg.read_output(outdir=str(out_dir))
        except (IOError, ValueError):
            continue
        X = np.ma.filled(fg.X, np.nan).reshape(-1, order="F")
        Y = np.ma.filled(fg.Y, np.nan).reshape(-1, order="F")
        B = np.ma.filled(fg.B, np.nan).reshape(-1, order="F")
        return {"X": X, "Y": Y, "B": B}
    return None


def plot_transects(mw_label: str, run_root: Path, geometry: Dict[str, np.ndarray], fgno: int) -> None:
    runs = list(iter_run_dirs(run_root))
    if not runs:
        return

    B0 = geometry["B"]
    Y = geometry["Y"]

    fig, ax = plt.subplots(figsize=(14, 5))

    for runno, run_dir in runs[:MAX_RUNS_PER_PLOT]:
        hmax_path = run_dir / f"hmax_fg{fgno}.npy"
        if not hmax_path.exists():
            print(f"[WARN] Missing {hmax_path}; skipping run {runno}")
            continue
        hmax = np.load(hmax_path)
        if hmax.shape[0] != B0.shape[0]:
            print(f"[WARN] Shape mismatch for {hmax_path}; skipping run {runno}")
            continue
        eta = np.where(hmax > 0.0, B0 + hmax, np.nan)
        ax.plot(Y, eta, color=ETA_COLOR, alpha=0.4)

    shoreline = np.where(B0 < 0.0, B0, np.nan)
    ax.fill_between(Y, shoreline, 0.0, color=(0.6, 0.6, 1.0), zorder=0)
    ax.plot(Y, B0, "g", linewidth=1.0, label="Topography")

    y_mid = 0.5 * (np.nanmin(Y) + np.nanmax(Y))
    ax.set_aspect(1.0 / math.cos(math.radians(y_mid)))
    ax.set_xlabel("Latitude (deg)")
    ax.set_ylabel("Elevation (m)")
    ax.grid(True)
    ax.set_title(f"Mw {mw_label} | fgmax {fgno} | first {min(len(runs), MAX_RUNS_PER_PLOT)} runs")

    plot_path = PNG_ROOT / f"Mw{mw_label.replace('.', '')}_fg{fgno}_transects.png"
    fig.savefig(plot_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Created {plot_path}")


def main() -> None:
    mw_entries = discover_mw_runs()
    if not mw_entries:
        print(f"[WARN] No Mw run directories found in {ALL_RUNS_DIR}")
        return

    for mw_label, run_root in mw_entries:
        raw_root = OUTPUT_ROOT / f"Mw_{mw_label.replace('.', '')}"
        for fgno in TRANSECT_FGMAX:
            geometry = load_fgmax_geometry(raw_root, fgno)
            if geometry is None:
                print(f"[WARN] Unable to load fgmax geometry for Mw {mw_label}, grid {fgno}")
                continue
            plot_transects(mw_label, run_root, geometry, fgno)


if __name__ == "__main__":
    main()
