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
import multiprocessing as mp
import os
import random
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
DATA_DIR = os.path.normpath(os.path.join(driver_home, '..', 'DataFiles'))
SCENARIO_DIR = os.path.join(driver_home, 'scenario')
DEFAULT_SCENARIO_FILE = os.path.join(SCENARIO_DIR, 'scenario_pts.txt')
OUTPUT_BASE = os.path.normpath(os.path.join(driver_home, '..', 'geoclaw_output'))

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
        default="test",
        help="Run random test subset ('test') or all scenarios ('all').",
    )
    parser.add_argument(
        "--n-test",
        type=int,
        default=5,
        help="Number of random cases to run in test mode.",
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
        help="Optional random seed for reproducible sampling in test mode.",
    )
    parser.add_argument(
        "--scenario-file",
        default=DEFAULT_SCENARIO_FILE,
        help="Path to the scenario points file (default: scenario/scenario_pts.txt).",
    )
    return parser.parse_args()


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

    if args.seed is not None:
        random.seed(args.seed)

    # ==========================================================================
    # Build base configuration once (for informational output / caching)
    # ==========================================================================
    get_fault()  # prints fault subdivision information

    scn_path = os.path.abspath(args.scenario_file)
    try:
        scn = np.loadtxt(scn_path)
    except OSError as err:
        print(f"Error reading scenario file {scn_path}: {err}")
        sys.exit(1)

    if isinstance(scn, np.ndarray):
        if scn.ndim == 0:
            scn_list = [float(scn)]
        elif scn.ndim == 1:
            scn_list = [scn.tolist()]
        else:
            scn_list = scn.tolist()
    else:
        scn_list = [scn]

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
