"""
Generate fgmax plots and summary statistics for all GeoClaw runs.

The script scans the `geoclaw_output` directory for Mw-specific runs,
produces contour plots of maximum tsunami depth for each run, and
records the maximum onshore runup and first arrival time.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from clawpack.geoclaw import fgmax_tools  # noqa: E402
from clawpack.geoclaw.geoplot import discrete_cmap_1  # noqa: E402
from clawpack.geoclaw.topotools import Topography  # noqa: E402


FGMAX_GRID_NO = 3
TARGET_MW_VALUES = ("8.0", "8.6", "8.8", "9.0")

BASE_DIR = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = BASE_DIR / "geoclaw_output"

TOPOGRAPHY_PATH = BASE_DIR / "DataFiles" / "CCtopo_1_3sec.tt3"


@dataclass
class RunSummary:
    """Container for per-run summary metrics."""

    magnitude: str
    run_name: str
    max_onshore_runup_m: Optional[float]
    max_onshore_depth_m: Optional[float]
    first_arrival_time_s: Optional[float]
    runup_lon: Optional[float]
    runup_lat: Optional[float]
    arrival_lon: Optional[float]
    arrival_lat: Optional[float]


def mw_directory_name(mw: str) -> str:
    """Translate Mw string (e.g. '8.6') to directory name."""
    code = mw.replace(".", "")
    return f"Mw_{code}"


def iter_run_directories() -> Iterable[Tuple[str, Path]]:
    """Yield (Mw_label, run_directory) pairs for configured magnitudes."""
    for mw_value in TARGET_MW_VALUES:
        mw_dir = OUTPUT_ROOT / mw_directory_name(mw_value)
        if not mw_dir.is_dir():
            continue
        for run_dir in sorted(mw_dir.glob("run_*")):
            if run_dir.is_dir():
                yield mw_value, run_dir


def load_fgmax_grid(run_dir: Path) -> Optional[fgmax_tools.FGmaxGrid]:
    """Load FGmax data for the requested run. Returns None if missing."""
    data_file = run_dir / "fgmax_grids.data"
    out_dir = run_dir / "_output"
    fgmax_file = out_dir / f"fgmax{FGMAX_GRID_NO:04d}.txt"

    if not data_file.exists() or not fgmax_file.exists():
        return None

    fg = fgmax_tools.FGmaxGrid()
    fg.read_fgmax_grids_data(FGMAX_GRID_NO, data_file=str(data_file))
    fg.read_output(outdir=str(out_dir))
    return fg


def compute_summary(fg: fgmax_tools.FGmaxGrid, mw_value: str, run_dir: Path) -> RunSummary:
    """Compute maximum onshore runup, depth, and first arrival time."""
    h = np.array(np.ma.filled(fg.h, np.nan))
    topo = np.array(np.ma.filled(fg.B, np.nan))
    eta = h + topo
    arrival = np.array(np.ma.filled(fg.arrival_time, np.nan))

    onshore_mask = topo > 0.0

    max_runup, runup_idx = _nanmax_with_index(eta, onshore_mask)
    max_depth, _ = _nanmax_with_index(h, onshore_mask)

    arrival_mask = onshore_mask & np.isfinite(arrival) & (arrival > 0.0)
    first_arrival, arrival_idx = _nanmin_with_index(arrival, arrival_mask)

    runup_lon, runup_lat = _coords_from_index(fg.X, fg.Y, runup_idx)
    arrival_lon, arrival_lat = _coords_from_index(fg.X, fg.Y, arrival_idx)

    return RunSummary(
        magnitude=mw_value,
        run_name=run_dir.name,
        max_onshore_runup_m=max_runup,
        max_onshore_depth_m=max_depth,
        first_arrival_time_s=first_arrival,
        runup_lon=runup_lon,
        runup_lat=runup_lat,
        arrival_lon=arrival_lon,
        arrival_lat=arrival_lat,
    )


def _nanmax_with_index(values: np.ndarray, mask: np.ndarray) -> Tuple[Optional[float], Optional[Tuple[int, ...]]]:
    """Return (max_value, index) considering mask; both None if no valid values."""
    masked = np.where(mask, values, np.nan)
    if np.all(np.isnan(masked)):
        return None, None
    flat_index = int(np.nanargmax(masked))
    index = np.unravel_index(flat_index, masked.shape)
    return float(masked[index]), index


def _nanmin_with_index(values: np.ndarray, mask: np.ndarray) -> Tuple[Optional[float], Optional[Tuple[int, ...]]]:
    """Return (min_value, index) considering mask; both None if no valid values."""
    masked = np.where(mask, values, np.nan)
    if np.all(np.isnan(masked)):
        return None, None
    flat_index = int(np.nanargmin(masked))
    index = np.unravel_index(flat_index, masked.shape)
    return float(masked[index]), index


def _coords_from_index(
    x: np.ndarray, y: np.ndarray, index: Optional[Tuple[int, ...]]
) -> Tuple[Optional[float], Optional[float]]:
    """Translate array indices to physical coordinates."""
    if index is None:
        return None, None
    return float(x[index]), float(y[index])


def plot_run(fg: fgmax_tools.FGmaxGrid, mw_value: str, run_dir: Path) -> None:
    """Create and save the fgmax plot for a single run."""
    h = np.ma.filled(fg.h, np.nan)

    clines = [0.01] + list(np.arange(0.5, 10.5, 0.5))
    colors = discrete_cmap_1(clines)

    fig, ax = plt.subplots(figsize=(10, 8))

    topo = _load_topography()
    contour_levels = np.arange(0.0, 21.0, 2.0)
    ax.contour(topo.X, topo.Y, topo.Z, levels=contour_levels, colors="#888888", linewidths=0.5)

    contourf = ax.contourf(fg.X, fg.Y, h, levels=clines, colors=colors, alpha=0.95, extend="max")
    cbar = fig.colorbar(contourf, ax=ax)
    cbar.set_label("Maximum water depth (m)")

    y0 = 0.5 * (np.nanmin(fg.Y) + np.nanmax(fg.Y))
    ax.set_aspect(1.0 / math.cos(math.radians(y0)))
    ax.ticklabel_format(useOffset=False)
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")

    ax.set_title(f"{run_dir.name} | Mw {mw_value}")

    x1 = float(np.nanmin(fg.X))
    x2 = float(np.nanmax(fg.X))
    y1 = float(np.nanmin(fg.Y))
    y2 = float(np.nanmax(fg.Y))
    ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], "g-")

    plot_dir = run_dir / "_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_path = plot_dir / f"{run_dir.name}_Mw{mw_value.replace('.', '')}_fgmax.png"
    fig.savefig(plot_path, bbox_inches="tight")
    plt.close(fig)


_cached_topography: Optional[Topography] = None


def _load_topography() -> Topography:
    """Load and cache the topography data used for context plots."""
    global _cached_topography  # noqa: PLW0603
    if _cached_topography is None:
        topo = Topography()
        topo.read(str(TOPOGRAPHY_PATH), 3)
        _cached_topography = topo
    return _cached_topography


def write_summary_csv(summaries: Iterable[RunSummary]) -> None:
    """Persist run summaries to CSV in the output root."""
    csv_path = OUTPUT_ROOT / "fgmax_run_summary.csv"
    fieldnames = [
        "Mw",
        "run",
        "max_onshore_runup_m",
        "max_onshore_depth_m",
        "first_arrival_time_s",
        "runup_lon",
        "runup_lat",
        "arrival_lon",
        "arrival_lat",
    ]
    with csv_path.open("w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for item in summaries:
            writer.writerow(
                {
                    "Mw": item.magnitude,
                    "run": item.run_name,
                    "max_onshore_runup_m": _format_optional_number(item.max_onshore_runup_m),
                    "max_onshore_depth_m": _format_optional_number(item.max_onshore_depth_m),
                    "first_arrival_time_s": _format_optional_number(item.first_arrival_time_s),
                    "runup_lon": _format_optional_number(item.runup_lon),
                    "runup_lat": _format_optional_number(item.runup_lat),
                    "arrival_lon": _format_optional_number(item.arrival_lon),
                    "arrival_lat": _format_optional_number(item.arrival_lat),
                }
            )


def _format_optional_number(value: Optional[float]) -> str:
    """Format optional float for CSV output."""
    if value is None or np.isnan(value):
        return ""
    return f"{value:.4f}"


def main() -> None:
    summaries = []
    for mw_value, run_dir in iter_run_directories():
        fg = load_fgmax_grid(run_dir)
        if fg is None:
            continue
        plot_run(fg, mw_value, run_dir)
        summary = compute_summary(fg, mw_value, run_dir)
        summaries.append(summary)
        runup_val = summary.max_onshore_runup_m if summary.max_onshore_runup_m is not None else float("nan")
        depth_val = summary.max_onshore_depth_m if summary.max_onshore_depth_m is not None else float("nan")
        arrival_val = summary.first_arrival_time_s if summary.first_arrival_time_s is not None else float("nan")
        print(
            f"{mw_value} {run_dir.name}: runup={runup_val:.3f} m, "
            f"depth={depth_val:.3f} m, arrival={arrival_val:.3f} s"
        )

    if summaries:
        write_summary_csv(summaries)
        print(f"Wrote summary CSV to {OUTPUT_ROOT / 'fgmax_run_summary.csv'}")


if __name__ == "__main__":
    main()
