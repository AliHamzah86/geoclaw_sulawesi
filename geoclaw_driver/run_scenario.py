"""
================================================

N O R T H   S U L A W E S I   T S U N A M I  
configuration : PTHA scenario for north sulawesi

================================================

run_PTHA for Sulawesi scenario:
    
"""

import argparse
from clawpack.geoclaw import dtopotools
import copy
import json
import multiprocessing as mp
import os
import random
import shutil
import sys
import time
import numpy as np

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

try:
    import rcrom
except:
    raise Exception("*** rcrom.py not in path: set PYTHONPATH")

# Scratch directory for storing topo and dtopo files:
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')
driver_home = os.getcwd()      # directory where all runs will be done

# Shared locations for input/output assets
DATA_DIR = os.path.join(driver_home, 'DataFiles')
SCENARIO_DIR = DATA_DIR
DEFAULT_SCENARIO_FILE = os.path.join(SCENARIO_DIR, 'scenario_pts.txt')
OUTPUT_BASE = os.path.normpath(os.path.join(driver_home, '..', 'geoclaw_output'))
DEFAULT_RANDOM_SEED = 12345
DEFAULT_COMPLETED_LOG = os.path.join(OUTPUT_BASE, 'completed_runs.log')
OUTPUT_KEEP_DIRS = {'common'}
DEFAULT_CONFIG_FILE = os.path.join(driver_home, 'run_config.json')

CONFIG_CACHE = {}

# Supported earthquake magnitudes and their output folders
MAGNITUDES = [8.0, 8.6, 8.8, 9.0]
MAG_OUTPUT_DIR = {
    8.0: 'Mw_80',
    8.6: 'Mw_86',
    8.8: 'Mw_88',
    9.0: 'Mw_90',
}


def magnitude_tag(mw):
    return MAG_OUTPUT_DIR[mw]


def format_run_tag(mw, scenario_index, run_id):
    mag_key = magnitude_tag(mw)
    return f"[{mag_key} run_{scenario_index} -> RUN-{run_id:03d}]"

# ==============================================================================
# setrun, setgeo for the grid runs are defined in setrun.py
#
#    setrun
#    setgeo
#
# these set as the default template, then the iteration function for the
# GeoClawInput class is used to appropriately change the settings,
# e.g., fine grid runs, earthquake magnitudes, run to final time, etc.
# ==============================================================================

from setrun import setrun, setgeo

# ==============================================================================
# Fault configuration constants and helpers
# ==============================================================================

LSTRIKE = 130e3
LDIP = 40e3
MAX_DEPTH = 20000.
PHI_PLATE = 60.0

COLUMN_MAP = {"longitude": 1, "latitude": 2, "depth": 3, "strike": 4,
              "length": 5, "width": 6, "dip": 7}
DEFAULTS = {'rake': 90, 'slip': 1.0}
COORDINATE_SPECIFICATION = 'top center'
INPUT_UNITS = {'slip': 'm', 'depth': 'km', 'length': 'km', 'width': 'km'}
RUPTURE_TYPE = 'static'
SKIPROWS = 1
DELIMITER = ','

BASE_FAULT = None


def depth_taper(d, max_depth=MAX_DEPTH):
    return 1. - np.exp((d - max_depth) * 5. / max_depth)


def build_fault():
    """Construct the subdivided Cascadia fault used for all runs."""
    fault = dtopotools.Fault()
    fault_path = os.path.join(DATA_DIR, 'CSZe01.csv')
    fault.read(fault_path, COLUMN_MAP, COORDINATE_SPECIFICATION,
               RUPTURE_TYPE, SKIPROWS, DELIMITER, INPUT_UNITS, DEFAULTS)
    print(f"There are {len(fault.subfaults)} subfaults")

    for subfault in fault.subfaults:
        subfault.longitude -= 360.0  # adjust to W coordinates

    # Restrict to southern portion (first 8 subfaults)
    fault.subfaults = fault.subfaults[:8]

    new_subfaults = []
    for subfault in fault.subfaults:
        subfault.rake = subfault.strike - PHI_PLATE - 180.0
        nstrike = int(subfault.length / 8000)
        ndip = int(subfault.width / 8000)
        subdivided = dtopotools.SubdividedPlaneFault(subfault, nstrike, ndip)
        new_subfaults.extend(subdivided.subfaults)

    new_fault = dtopotools.Fault(subfaults=new_subfaults)
    print(f"Subdivided fault has {len(new_fault.subfaults)} subfaults")
    return new_fault


def get_fault():
    """Return the cached base fault geometry, constructing it if needed."""
    global BASE_FAULT
    if BASE_FAULT is None:
        BASE_FAULT = build_fault()
    return BASE_FAULT


def format_duration(seconds):
    """Return a compact human-readable string for a duration in seconds."""
    seconds = float(seconds)
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = seconds / 60.0
    if minutes < 60:
        return f"{minutes:.1f}m"
    hours = minutes / 60.0
    return f"{hours:.1f}h"


def describe_run(run_id, M):
    """Return metadata (grid, Mw, scenario index, rundir) for a given run_id."""
    total_runs = len(MAGNITUDES) * M
    metadata = {
        "run_id": run_id,
        "mw": 0.0,
        "scenario": None,
        "rundir": None,
        "label": "",
    }

    if run_id < 0 or run_id >= total_runs:
        raise ValueError(f"run_id {run_id} exceeds expected range (0-{total_runs - 1})")

    scenario_index = run_id % M
    mag_index = run_id // M
    Mw = MAGNITUDES[mag_index]

    metadata["scenario"] = scenario_index
    metadata["mw"] = Mw
    metadata["mag_tag"] = magnitude_tag(Mw)
    metadata["tag"] = format_run_tag(Mw, scenario_index, run_id)

    label = f"Mw {Mw:.1f}"
    out_dir = os.path.join(OUTPUT_BASE, metadata["mag_tag"])
    metadata["label"] = label
    metadata["rundir"] = os.path.join(out_dir, f"run_{scenario_index}")

    return metadata


def parse_cli_args():
    parser = argparse.ArgumentParser(
        description="Run GeoClaw scenarios in test or full mode."
    )
    parser.add_argument(
        "--mode",
        choices=["test", "all"],
        default=None,
        help="Run random test subset ('test') or all scenarios ('all'). (default: configuration file)",
    )
    parser.add_argument(
        "--n-test",
        type=int,
        default=None,
        help="Number of random cases to run in test mode. (default: configuration file)",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=None,
        help="Maximum number of parallel worker processes (default: min(8, CPU count)).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible sampling in test mode. (default: configuration file)",
    )
    parser.add_argument(
        "--scenario-file",
        default=None,
        help="Path to the scenario points file. (default: configuration file)",
    )
    parser.add_argument(
        "--resume",
        dest="resume",
        action="store_true",
        default=None,
        help="Skip run IDs recorded in the completed log and append new successes.",
    )
    parser.add_argument(
        "--no-resume",
        dest="resume",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--completed-log",
        default=None,
        help="Path to the completed-runs log. (default: configuration file)",
    )
    parser.add_argument(
        "--max-cases",
        type=int,
        default=None,
        help="Limit the number of runs executed this invocation (after filtering).",
    )
    parser.add_argument(
        "--clean-output",
        dest="clean_output",
        action="store_true",
        default=None,
        help="Remove all contents of geoclaw_output except the 'common' directory before running.",
    )
    parser.add_argument(
        "--no-clean-output",
        dest="clean_output",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--run-after-clean",
        action="store_true",
        help="If provided with --clean-output, continue to run scenarios after cleaning.",
    )
    parser.add_argument(
        "--show-options",
        action="store_true",
        help="Display available command-line options and exit.",
    )
    parser.add_argument(
        "--status-only",
        action="store_true",
        help="Report completion status using the scenario file and completed log, then exit.",
    )
    parser.add_argument(
        "--config-file",
        default=DEFAULT_CONFIG_FILE,
        help=f"Path to the configuration JSON file (default: {DEFAULT_CONFIG_FILE}).",
    )
    parser.add_argument(
        "--write-config",
        action="store_true",
        help="Write configuration file with current settings before running.",
    )
    parser.add_argument(
        "--info-config",
        action="store_true",
        help="Display the resolved configuration (including overrides) and exit.",
    )
    args = parser.parse_args()
    setattr(args, "_parser", parser)
    return args


def chunk_run_ids(run_ids, max_chunks):
    """Split run_ids into at most max_chunks contiguous groups."""
    if max_chunks <= 0:
        return []

    run_ids_sorted = sorted(run_ids)
    if not run_ids_sorted:
        return []

    max_chunks = min(max_chunks, len(run_ids_sorted))
    splits = np.array_split(run_ids_sorted, max_chunks)
    return [list(map(int, chunk)) for chunk in splits if len(chunk) > 0]


def select_balanced_test_runs(M, ntest):
    """Select run IDs for test mode balancing across magnitudes."""
    num_mags = len(MAGNITUDES)
    ntest = max(1, ntest)
    counts = [0] * num_mags
    base = ntest // num_mags
    remainder = ntest % num_mags

    for idx in range(num_mags):
        counts[idx] = base

    if remainder > 0:
        mag_indices = list(range(num_mags))
        random.shuffle(mag_indices)
        for idx in mag_indices[:remainder]:
            counts[idx] += 1

    selected = []
    per_mag_counts = {mw: 0 for mw in MAGNITUDES}
    for mag_index, count in enumerate(counts):
        mag_ids = list(range(mag_index * M, (mag_index + 1) * M))
        if not mag_ids or count <= 0:
            continue
        actual = min(count, len(mag_ids))
        sampled = random.sample(mag_ids, actual)
        selected.extend(sampled)
        per_mag_counts[MAGNITUDES[mag_index]] = len(sampled)

    return sorted(selected[:ntest]), per_mag_counts


def load_completed_run_ids(log_path):
    """Return a set of completed run IDs recorded in log_path."""
    try:
        with open(log_path, 'r', encoding='ascii') as handle:
            return {int(line.strip()) for line in handle if line.strip()}
    except FileNotFoundError:
        return set()


def append_completed_run_ids(log_path, completed_ids):
    """Append completed run IDs to the log, creating it if necessary."""
    if not completed_ids:
        return
    directory = os.path.dirname(log_path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    # Convert to sorted unique entries before appending for readability.
    unique_ids = sorted(set(completed_ids))
    with open(log_path, 'a', encoding='ascii') as handle:
        for run_id in unique_ids:
            handle.write(f"{run_id}\n")


def default_config():
    """Return the default configuration dictionary."""
    return {
        "mode": "test",
        "n_test": 5,
        "seed": DEFAULT_RANDOM_SEED,
        "processes": None,
        "max_cases": None,
        "resume": False,
        "clean_output": False,
        "scenario_file": DEFAULT_SCENARIO_FILE,
        "completed_log": DEFAULT_COMPLETED_LOG,
        "output_base": OUTPUT_BASE,
        "data_dir": DATA_DIR,
        "topography": [
            {
                "topotype": 3,
                "minlevel": 1,
                "maxlevel": 4,
                "t1": 0.0,
                "t2": 1.0e10,
                "path": os.path.join(DATA_DIR, 'etopo1_-130_-124_38_45_1min.asc'),
            },
        ],
        "fgmax_files": [
            {"grid_id": 1, "path": os.path.join(DATA_DIR, 'fgmax1_coarse.txt')},
            {"grid_id": 2, "path": os.path.join(DATA_DIR, 'fgmax2_coarse.txt')},
            {"grid_id": 3, "path": os.path.join(DATA_DIR, 'fgmax3_coarse.txt')},
        ],
    }


def load_config_file(path):
    """Load configuration from JSON file, returning {} if not found."""
    try:
        with open(path, 'r', encoding='utf-8') as handle:
            return json.load(handle)
    except FileNotFoundError:
        return {}
    except json.JSONDecodeError as err:
        raise RuntimeError(f"Error parsing configuration file {path}: {err}")


def write_config_file(path, config):
    """Write configuration dictionary to JSON."""
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(config, handle, indent=2, sort_keys=True)
        handle.write("\n")


def merge_config(base_config, override):
    """Return a shallow merge of base_config with override dict."""
    merged = dict(base_config)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_config(merged.get(key, {}), value)
        else:
            merged[key] = value
    return merged


def _normalize_path(value):
    if value is None:
        return None
    return os.path.abspath(os.path.expanduser(value))


def normalize_config_paths(config):
    """Return configuration with key paths converted to absolute form."""
    normalized = dict(config)
    for key in ["scenario_file", "completed_log", "output_base", "data_dir"]:
        if key in normalized and normalized[key] is not None:
            normalized[key] = _normalize_path(normalized[key])

    if "topography" in normalized:
        topo_entries = []
        for entry in normalized["topography"]:
            topo_copy = dict(entry)
            if "path" in topo_copy and topo_copy["path"] is not None:
                topo_copy["path"] = _normalize_path(topo_copy["path"])
            topo_entries.append(topo_copy)
        normalized["topography"] = topo_entries

    if "fgmax_files" in normalized:
        fgmax_entries = []
        for entry in normalized["fgmax_files"]:
            fgmax_copy = dict(entry)
            if "path" in fgmax_copy and fgmax_copy["path"] is not None:
                fgmax_copy["path"] = _normalize_path(fgmax_copy["path"])
            fgmax_entries.append(fgmax_copy)
        normalized["fgmax_files"] = fgmax_entries

    return normalized


def ensure_config_paths(config):
    """Fallback to default locations if configured paths are missing."""
    updated = dict(config)

    data_dir = updated.get("data_dir")
    if data_dir and not os.path.isdir(data_dir):
        default_dir = os.path.join(driver_home, 'DataFiles')
        if os.path.isdir(default_dir):
            print(f"Config | data_dir {data_dir} not found, using {default_dir}")
            updated["data_dir"] = default_dir

    scenario_file = updated.get("scenario_file")
    if scenario_file and not os.path.exists(scenario_file):
        fallback = os.path.join(updated.get("data_dir", os.path.join(driver_home, 'DataFiles')), 'scenario_pts.txt')
        if os.path.exists(fallback):
            print(f"Config | scenario_file {scenario_file} not found, using {fallback}")
            updated["scenario_file"] = fallback

    completed_log = updated.get("completed_log")
    if completed_log:
        log_dir = os.path.dirname(completed_log)
        if log_dir and not os.path.isdir(log_dir):
            try:
                os.makedirs(log_dir, exist_ok=True)
            except OSError:
                pass

    topo_entries = updated.get("topography")
    if topo_entries:
        fixed_entries = []
        for entry in topo_entries:
            entry_copy = dict(entry)
            path = entry_copy.get("path")
            if path and not os.path.exists(path):
                alt_path = os.path.join(updated.get("data_dir", DATA_DIR), os.path.basename(path))
                if os.path.exists(alt_path):
                    print(f"Config | topography file {path} not found, using {alt_path}")
                    entry_copy["path"] = alt_path
            fixed_entries.append(entry_copy)
        updated["topography"] = fixed_entries

    fgmax_entries = updated.get("fgmax_files")
    if fgmax_entries:
        fixed_fgmax = []
        for entry in fgmax_entries:
            entry_copy = dict(entry)
            path = entry_copy.get("path")
            if path and not os.path.exists(path):
                basename = os.path.basename(path)
                candidates = [
                    os.path.join(updated.get("data_dir", DATA_DIR), basename),
                    os.path.join(driver_home, basename),
                ]
                for alt_path in candidates:
                    if os.path.exists(alt_path):
                        print(f"Config | FGmax file {path} not found, using {alt_path}")
                        entry_copy["path"] = alt_path
                        break
            fixed_fgmax.append(entry_copy)
        updated["fgmax_files"] = fixed_fgmax

    return updated


def config_overrides_from_cli(args):
    """Collect configuration overrides from provided CLI arguments."""
    overrides = {}
    if args.mode is not None:
        overrides["mode"] = args.mode
    if args.n_test is not None:
        overrides["n_test"] = args.n_test
    if args.seed is not None:
        overrides["seed"] = args.seed
    if args.processes is not None:
        overrides["processes"] = args.processes
    if args.max_cases is not None:
        overrides["max_cases"] = args.max_cases
    if args.resume is not None:
        overrides["resume"] = args.resume
    if args.clean_output is not None:
        overrides["clean_output"] = args.clean_output
    if args.scenario_file is not None:
        overrides["scenario_file"] = args.scenario_file
    if args.completed_log is not None:
        overrides["completed_log"] = args.completed_log
    return overrides


def derive_runtime_settings(args, config):
    """Blend CLI arguments with configuration to produce effective settings."""
    settings = {}
    mode_value = config.get("mode", "test")
    settings["mode"] = args.mode if args.mode is not None else (mode_value if mode_value in {"test", "all"} else "test")

    ntest_value = config.get("n_test", 5)
    if args.n_test is not None:
        settings["n_test"] = args.n_test
    else:
        settings["n_test"] = ntest_value if ntest_value is not None else 5

    seed_value = config.get("seed", DEFAULT_RANDOM_SEED)
    settings["seed"] = args.seed if args.seed is not None else seed_value
    settings["processes"] = args.processes if args.processes is not None else config.get("processes")
    settings["max_cases"] = args.max_cases if args.max_cases is not None else config.get("max_cases")
    settings["resume"] = args.resume if args.resume is not None else bool(config.get("resume", False))
    settings["clean_output"] = args.clean_output if args.clean_output is not None else bool(config.get("clean_output", False))
    scenario_value = config.get("scenario_file", DEFAULT_SCENARIO_FILE)
    if scenario_value is None:
        scenario_value = DEFAULT_SCENARIO_FILE
    settings["scenario_file"] = _normalize_path(
        args.scenario_file if args.scenario_file is not None else scenario_value
    )

    completed_value = config.get("completed_log")
    if completed_value is None:
        completed_value = os.path.join(config.get("output_base", OUTPUT_BASE), 'completed_runs.log')
    settings["completed_log"] = _normalize_path(
        args.completed_log if args.completed_log is not None else completed_value
    )
    settings["config_file"] = _normalize_path(args.config_file)
    return settings


def apply_runtime_paths_from_config(config):
    """Update module-level directory paths from configuration."""
    global DATA_DIR, OUTPUT_BASE
    if config.get("data_dir"):
        DATA_DIR = config["data_dir"]
    if config.get("output_base"):
        OUTPUT_BASE = config["output_base"]


def clean_output_directory(base_dir):
    """Remove everything in base_dir except directories listed in OUTPUT_KEEP_DIRS."""
    removed = []
    if not os.path.isdir(base_dir):
        return removed
    for entry in os.listdir(base_dir):
        if entry in OUTPUT_KEEP_DIRS:
            continue
        path = os.path.join(base_dir, entry)
        try:
            if os.path.isdir(path) and not os.path.islink(path):
                shutil.rmtree(path)
                removed.append(f"{entry}/")
            else:
                os.remove(path)
                removed.append(entry)
        except FileNotFoundError:
            continue
    return removed


def report_run_status(M, known_completed, flat_results):
    """Print concise per-case completion info and any encountered errors."""
    print("\nRun Status")
    for mag_index, Mw in enumerate(MAGNITUDES):
        tag = magnitude_tag(Mw)
        run_ids = range(mag_index * M, (mag_index + 1) * M)
        completed = sorted(rid % M for rid in run_ids if rid in known_completed)
        remaining = sorted(rid % M for rid in run_ids if rid not in known_completed)
        completed_text = ", ".join(f"run_{idx}" for idx in completed) if completed else "none"
        print(f"  [{tag}] completed: {completed_text}")
        print(f"  [{tag}] remaining: {len(remaining)} case(s)")

    failures = [(run_id, status) for run_id, status, _ in flat_results if status.startswith("FAIL")]
    if failures:
        short_failures = ", ".join(f"RUN-{run_id:03d}" for run_id, _ in failures)
        print(f"  Errors: {short_failures}")
    else:
        print("  Errors: none")


def load_scenario_points(path):
    """Return list of scenario points parsed from a file."""
    try:
        scn = np.loadtxt(path)
    except OSError as err:
        raise RuntimeError(f"Error reading scenario file {path}: {err}")

    if isinstance(scn, np.ndarray):
        if scn.ndim == 0:
            return [float(scn)]
        if scn.ndim == 1:
            return [scn.tolist()]
        return scn.tolist()
    return [scn]


# ==============================================================================
# Iterator for the runs 
#
#       iter_fun() is the iterator function for the GeoClawInput class 
#       run parameters (grid-size, earthquake parameters) as well as 
#       run_ids, rundirs, etc. are changed here
#       
# ==============================================================================

def iter_fun(self):
    r"""
    This function will be used to iterate GeoClawInput
    Total runs: len(MAGNITUDES) * len(scenarios)
    """

    # helper functions to set slip distribution
    from clawpack.geoclaw import fgmax_tools


    run_id = self._run_id
    etopo_dir = DATA_DIR
    topodir = DATA_DIR

    # load input info
    if self._input_info == None:
        scn_fname = DEFAULT_SCENARIO_FILE
        scn = np.loadtxt(scn_fname)
        scn_list = scn.tolist()
    else:
        scn_list = self._input_info
    
    # total number of runs
    M = len(scn_list)
    N = len(MAGNITUDES) * M

    if run_id == N:
        raise StopIteration()

    else:
        run_id_mod = run_id % M
        mag_index = run_id // M
        self.KL_Mw_desired = MAGNITUDES[mag_index]
        run_tag = format_run_tag(self.KL_Mw_desired, run_id_mod, run_id)

        #=========================
        # set grid configurations
        #
        t_shelf = 0.   # time approaching continental slope
        t_harbor = 0.  # time approaching harbor

        self._rundata = setrun(setgeo)   # includes fgmax setup

        self._rundata.amrdata.amr_levels_max = 4
        # grid run = 10sec
        # dx = 30min, 5min, 1min, 10sec
        self._rundata.amrdata.refinement_ratios_x = [6, 5, 6]
        self._rundata.amrdata.refinement_ratios_y = [6, 5, 6]
        self._rundata.amrdata.refinement_ratios_t = [6, 5, 6]

        # add topography
        topofiles = self._rundata.topo_data.topofiles
        # for topography, append lines of the form
        #    [topotype, minlevel, maxlevel, t1, t2, fname]
        topofiles[:] = []  # clear existing entries but keep reference

        topo_config = CONFIG_CACHE.get("topography", []) if CONFIG_CACHE else []
        if topo_config:
            for entry in topo_config:
                if "path" not in entry:
                    raise ValueError("Topography configuration entries must include a 'path'.")
                topo_record = [
                    entry.get("topotype", 3),
                    entry.get("minlevel", 1),
                    entry.get("maxlevel", 4),
                    entry.get("t1", 0.0),
                    entry.get("t2", 1.0e10),
                    entry["path"],
                ]
                topofiles.append(topo_record)
        else:
            topofiles.append([3, 1, 4, 0., 1.e10,
                    os.path.join(etopo_dir, 'etopo1_-130_-124_38_45_1min.asc')])
            topofiles.append([3, 3, 4, 0., 1.e10,
                    os.path.join(topodir, 'cc-1sec.asc')])

        # add regions
        regions = self._rundata.regiondata.regions 
        # between shelf and CC 
        regions = []
        regions.append(
                [2, 3, t_shelf, 1e9, -125, -124.05, 40.5, 43]) 
        regions.append(
                [3, 4, t_harbor, 1e9, -124.26, -124.14, 41.67,   41.79])
        regions.append(
                [4, 4, t_harbor, 1e9, -124.218,-124.17, 41.7345, 41.77])
        self._rundata.regiondata.regions = regions   # ✅ put it back

        # check that topo files exist
        for _, _, _, _, _, f in topofiles:
            if not os.path.exists(f):
                print(f"{run_tag} | missing topography file: {f}")
                raise FileNotFoundError(f)

        fgmax_config = CONFIG_CACHE.get("fgmax_files", []) if CONFIG_CACHE else []
        if fgmax_config:
            fgmax_grids = self._rundata.fgmax_data.fgmax_grids
            fgmax_grids[:] = []
            for entry in fgmax_config:
                if "path" not in entry or "grid_id" not in entry:
                    raise ValueError("FGmax configuration entries must include 'grid_id' and 'path'.")
                grid = fgmax_tools.FGmaxGrid()
                if not os.path.exists(entry["path"]):
                    print(f"{run_tag} | missing FGmax file: {entry['path']}")
                    raise FileNotFoundError(entry["path"])
                grid.read_fgmax_grids_data(entry["grid_id"], entry["path"])
                fgmax_grids.append(grid)

        #
        # set slip distribution
        #
        # run_id_mod = run_id - 100*(run_id/100)
        m = scn_list[run_id_mod]
        self.set_KL_slip(m)
    
        dir_grid_Mw = os.path.join(OUTPUT_BASE, magnitude_tag(self.KL_Mw_desired))
        os.makedirs(dir_grid_Mw, exist_ok=True)
        self._rundir = os.path.join(dir_grid_Mw, 'run_' + str(run_id_mod))
        
        # --- Compact progress info (after Mw assigned) ---
        if getattr(self, "_progress_enabled", True):
            print(f"{run_tag} scheduled ({run_id + 1}/{N})")
            
        self._run_id += 1
        
        return self


def create_configured_drom(scenario_points):
    """Create a Drom instance with GeoClaw input configured for this scenario."""
    drom = rcrom.Drom()
    drom.GeoClawInput.fault = copy.deepcopy(get_fault())
    drom.GeoClawInput.set_iter(iter_fun)
    drom.GeoClawInput.set_rundata(setrun=setrun, setgeo=setgeo)
    drom.GeoClawInput.KL_expand(Lstrike=LSTRIKE, Ldip=LDIP,
                distribution='Lognormal', tau=depth_taper,
                nterms=20, KL_Mw_desired=9.0)
    drom.GeoClawInput._input_info = scenario_points
    return drom


def run_cases_sequential(run_ids, driver_path, scenario_points):
    """Run the GeoClaw cases identified by run_ids sequentially."""
    run_ids = sorted(run_ids)
    if not run_ids:
        return []

    pending = set(run_ids)
    results = []

    try:
        os.chdir(driver_path)
    except OSError:
        pass

    drom = create_configured_drom(scenario_points)
    scenario_count = len(scenario_points) if scenario_points is not None else 0
    drom.GeoClawInput._progress_enabled = False

    for geoinput in drom.GeoClawInput:
        current_id = geoinput._run_id - 1

        if current_id not in pending:
            continue

        run_dir = geoinput._rundir
        Mw = geoinput.KL_Mw_desired
        run_id_mod = current_id % scenario_count if scenario_count else current_id
        run_tag = format_run_tag(Mw, run_id_mod, current_id)

        orig_geo_id = drom.GeoClawInput._run_id
        orig_drom_id = drom._run_id
        drom.GeoClawInput._run_id = current_id
        drom._run_id = current_id

        start = time.perf_counter()
        try:
            drom.evaluate_hdm(run_tag=run_tag)
            duration = time.perf_counter() - start
            print(f"{run_tag} | completed in {format_duration(duration)}")
            results.append((current_id, "OK", duration))
        except Exception as exc:
            duration = time.perf_counter() - start
            print(f"{run_tag} | failed: {exc} ({format_duration(duration)})")
            results.append((current_id, f"FAIL: {exc}", duration))
        finally:
            drom.GeoClawInput._run_id = orig_geo_id
            drom._run_id = orig_drom_id

        pending.remove(current_id)
        if not pending:
            break

    for remaining in sorted(pending):
        results.append((remaining, "SKIPPED", 0.0))

    return results


def run_case_worker(args):
    """Worker wrapper to run a set of cases with multiprocessing."""
    run_ids, driver_path, scenario_points = args
    return run_cases_sequential(run_ids, driver_path, scenario_points)


if __name__=='__main__':
    
    args = parse_cli_args()
    parser = getattr(args, "_parser", None)

    if getattr(args, "show_options", False):
        if parser is not None:
            print(parser.format_help())
        else:
            print("No option information available.")
        sys.exit(0)

    config_path = _normalize_path(getattr(args, "config_file", DEFAULT_CONFIG_FILE))
    args.config_file = config_path

    base_config = default_config()
    try:
        file_config = load_config_file(config_path)
    except RuntimeError as err:
        print(err)
        sys.exit(1)

    config_exists = os.path.exists(config_path)
    config = merge_config(base_config, file_config)
    config = normalize_config_paths(config)
    config = ensure_config_paths(config)

    if getattr(args, "write_config", False):
        overrides = config_overrides_from_cli(args)
        if overrides:
            config = merge_config(config, overrides)
            config = normalize_config_paths(config)
            config = ensure_config_paths(config)
        try:
            write_config_file(config_path, config)
            print(f"Config | wrote configuration to {config_path}")
        except OSError as err:
            print(f"Config | failed to write configuration {config_path}: {err}")
            sys.exit(1)
    elif not config_exists:
        try:
            write_config_file(config_path, config)
            print(f"Config | created default configuration at {config_path}")
        except OSError as err:
            print(f"Config | failed to create default configuration {config_path}: {err}")
            sys.exit(1)

    apply_runtime_paths_from_config(config)
    CONFIG_CACHE = config

    runtime_settings = derive_runtime_settings(args, config)
    for key, value in runtime_settings.items():
        setattr(args, key, value)

    DEFAULT_SCENARIO_FILE = args.scenario_file
    DEFAULT_COMPLETED_LOG = args.completed_log

    CONFIG_CACHE["scenario_file"] = args.scenario_file
    CONFIG_CACHE["completed_log"] = args.completed_log
    CONFIG_CACHE["data_dir"] = DATA_DIR
    CONFIG_CACHE["output_base"] = OUTPUT_BASE

    if getattr(args, "info_config", False):
        print("Configuration Summary")
        print(f"  config file    : {args.config_file}")
        print(f"  mode (default) : {CONFIG_CACHE.get('mode')}")
        print(f"  n_test (default): {CONFIG_CACHE.get('n_test')}")
        print(f"  seed (default) : {CONFIG_CACHE.get('seed')}")
        print(f"  data dir       : {CONFIG_CACHE.get('data_dir')}")
        print(f"  scenario file  : {CONFIG_CACHE.get('scenario_file')}")
        print(f"  completed log  : {CONFIG_CACHE.get('completed_log')}")
        print(f"  output base    : {CONFIG_CACHE.get('output_base')}")
        topography = CONFIG_CACHE.get("topography", [])
        if topography:
            print("  topography files:")
            for entry in topography:
                print(f"    - level {entry.get('minlevel')}–{entry.get('maxlevel')} ({entry.get('topotype')}): {entry.get('path')}")
        else:
            print("  topography files: none configured")
        fgmax_entries = CONFIG_CACHE.get("fgmax_files", [])
        if fgmax_entries:
            print("  FGmax grids:")
            for entry in fgmax_entries:
                print(f"    - grid {entry.get('grid_id')}: {entry.get('path')}")
        else:
            print("  FGmax grids: none configured")
        print("  runtime overrides this run:")
        for key in ["mode", "n_test", "seed", "processes", "max_cases", "resume", "clean_output"]:
            print(f"    {key}: {runtime_settings.get(key)}")
        sys.exit(0)

    print(f"Info | config file: {args.config_file}")
    print(f"Info | data dir: {DATA_DIR}")
    print(f"Info | scenario file: {args.scenario_file}")
    print(f"Info | output base: {OUTPUT_BASE}")

    if getattr(args, "status_only", False):
        try:
            scn_list = load_scenario_points(args.scenario_file)
        except RuntimeError as err:
            print(err)
            sys.exit(1)

        M = len(scn_list)
        if M == 0:
            print("Status only | scenario list is empty.")
            sys.exit(0)

        total_runs = len(MAGNITUDES) * M
        completed_run_ids = load_completed_run_ids(args.completed_log)
        print(f"Status only | total configured runs: {total_runs}")
        report_run_status(M, completed_run_ids, [])
        sys.exit(0)

    if args.clean_output:
        removed_entries = clean_output_directory(OUTPUT_BASE)
        if removed_entries:
            print(f"Clean output | removed {len(removed_entries)} item(s): {', '.join(removed_entries)}")
        else:
            print("Clean output | nothing to remove (directories already clean).")
        if not getattr(args, "run_after_clean", False):
            print("Clean output | skipping scenario execution.")
            sys.exit(0)

    if args.seed is not None:
        random.seed(args.seed)

    completed_run_ids = set()
    if args.resume:
        completed_run_ids = load_completed_run_ids(args.completed_log)
        if completed_run_ids:
            print(f"Resume enabled | loaded {len(completed_run_ids)} completed run(s).")
        else:
            print(f"Resume enabled | no completed runs recorded yet.")

    # ==========================================================================
    # Build base configuration once (for informational output / caching)
    # ==========================================================================
    get_fault()  # prints fault subdivision information

    scn_path = os.path.abspath(args.scenario_file)
    try:
        scn_list = load_scenario_points(scn_path)
    except RuntimeError as err:
        print(err)
        sys.exit(1)

    M = len(scn_list)
    if M == 0:
        print("Error: scenario list is empty.")
        sys.exit(1)

    total_runs = len(MAGNITUDES) * M
    all_run_ids = list(range(total_runs))

    if args.mode == "test":
        ntest = max(1, min(args.n_test, total_runs))
        selected_run_ids, mag_counts = select_balanced_test_runs(M, ntest)
        count_text = ", ".join(f"{magnitude_tag(mw)}: {cnt}" for mw, cnt in sorted(mag_counts.items()))
        print(f"Mode test | {len(selected_run_ids)} run(s) selected [{count_text}]")
    else:
        selected_run_ids = all_run_ids
        print(f"Mode all | scheduling {total_runs} run(s).")

    if args.resume and selected_run_ids:
        before = len(selected_run_ids)
        selected_run_ids = [rid for rid in selected_run_ids if rid not in completed_run_ids]
        skipped = before - len(selected_run_ids)
        if skipped:
            print(f"Resume filter | skipping {skipped} previously completed run(s).")
        elif completed_run_ids:
            print("Resume filter | all selected runs are new.")

    if args.max_cases is not None:
        if args.max_cases <= 0:
            print("Error: max-cases must be positive.")
            sys.exit(1)
        if len(selected_run_ids) > args.max_cases:
            print(f"Batch limit | reducing to first {args.max_cases} run(s).")
            selected_run_ids = selected_run_ids[:args.max_cases]
        else:
            print(f"Batch limit | {len(selected_run_ids)} run(s) within specified maximum.")

    required_dirs = {os.path.join(OUTPUT_BASE, 'common')}
    required_dirs.update(os.path.join(OUTPUT_BASE, name) for name in MAG_OUTPUT_DIR.values())
    for path in sorted(required_dirs):
        os.makedirs(path, exist_ok=True)

    if not selected_run_ids:
        print("Error: no runs selected.")
        sys.exit(1)

    print("Scheduled runs:")
    for rid in selected_run_ids:
        meta = describe_run(rid, M)
        print(f"  {meta['tag']} | directory {meta['rundir']}")

    cpu_count = mp.cpu_count()
    max_workers = args.processes if args.processes is not None else min(8, cpu_count)
    if max_workers <= 0:
        print("Error: number of processes must be positive.")
        sys.exit(1)

    run_chunks = chunk_run_ids(selected_run_ids, max_workers)
    if not run_chunks:
        print("Error: failed to partition run IDs for execution.")
        sys.exit(1)

    actual_workers = len(run_chunks)
    if actual_workers == 1:
        print(f"Execution mode: sequential (CPU count {cpu_count}).")
    else:
        print(f"Execution mode: {actual_workers} workers (CPU count {cpu_count}).")

    worker_args = [(chunk, driver_home, scn_list) for chunk in run_chunks]

    wall_start = time.perf_counter()
    if actual_workers == 1:
        flat_results = run_cases_sequential(worker_args[0][0], driver_home, scn_list)
    else:
        with mp.Pool(processes=actual_workers) as pool:
            worker_results = pool.map(run_case_worker, worker_args)
        flat_results = [item for sublist in worker_results for item in sublist]
    wall_total = time.perf_counter() - wall_start

    if args.resume:
        successful_ids = [run_id for run_id, status, _ in flat_results if status == "OK"]
        new_completed = sorted(set(successful_ids) - completed_run_ids)
        append_completed_run_ids(args.completed_log, new_completed)
        completed_run_ids.update(new_completed)
        if new_completed:
            print(f"Resume log | recorded {len(new_completed)} new run(s).")
        elif successful_ids:
            print("Resume log | all successful runs were already recorded.")

    ok = sum(1 for _, status, _ in flat_results if status == "OK")
    failures = [(run_id, status, duration) for run_id, status, duration in flat_results if status.startswith("FAIL")]
    skipped = [(run_id, status, duration) for run_id, status, duration in flat_results if status == "SKIPPED"]
    total_run_time = sum(duration for _, _, duration in flat_results)
    total_success_time = sum(duration for _, status, duration in flat_results if status == "OK")

    mag_total = {mw: 0 for mw in MAGNITUDES}
    mag_success = {mw: 0 for mw in MAGNITUDES}
    for run_id, status, _ in flat_results:
        meta = describe_run(run_id, M)
        mw = meta["mw"]
        mag_total[mw] += 1
        if status == "OK":
            mag_success[mw] += 1

    print("\nSummary")
    print(f"  Successful runs: {ok}")
    if failures:
        for run_id, status, duration in sorted(failures):
            meta = describe_run(run_id, M)
            print(f"  Failed {meta['tag']}: {status} ({format_duration(duration)})")
    else:
        print("  Failed runs: 0")
    if skipped:
        skipped_tags = [describe_run(run_id, M)['tag'] for run_id, _, _ in skipped]
        print(f"  Skipped runs: {', '.join(skipped_tags)}")
    print(f"  Total wall time: {format_duration(wall_total)}")
    print(f"  Sum of run durations: {format_duration(total_run_time)}")
    if ok:
        avg_duration = total_success_time / ok
        print(f"  Mean successful run: {format_duration(avg_duration)}")
    print("  Cases per magnitude:")
    for mw in MAGNITUDES:
        tag = magnitude_tag(mw)
        total = mag_total.get(mw, 0)
        success = mag_success.get(mw, 0)
        plural = "case" if total == 1 else "cases"
        print(f"    [{tag}]: {total} {plural} ({success} successful)")
    print("--------------------------------------------------")

    if args.resume:
        known_completed = set(completed_run_ids)
    else:
        known_completed = {run_id for run_id, status, _ in flat_results if status == "OK"}
    report_run_status(M, known_completed, flat_results)
