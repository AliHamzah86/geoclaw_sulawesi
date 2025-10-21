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
data_dir = os.path.join(driver_home, 'DataFiles')

# ==============================================================================
# setrun, setgeo for the coarse grid runs are defined in setrun.py
#
#    setrun_coarse
#    setgeo_coarse
#
# these set as the default template, then the iteration function for the
# GeoClawInput class is used to appropriately change the settings,
# e.g., fine grid runs, earthquake magnitudes, run to final time, etc.
# ==============================================================================

from setrun import setrun_coarse, setgeo_coarse

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

    fault.read(os.path.join(data_dir,'CSZe01.csv'), COLUMN_MAP, COORDINATE_SPECIFICATION,
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
    mw_cycle = [8.6, 8.8, 9.0, 9.2]
    total_non_b0 = 8 * M
    metadata = {
        "run_id": run_id,
        "grid": None,
        "mw": 0.0,
        "scenario": None,
        "rundir": None,
        "label": "",
    }

    if run_id < 0:
        raise ValueError("run_id must be non-negative")

    if run_id < 4 * M:
        metadata["grid"] = "coarse"
        metadata["scenario"] = run_id % M
        band = run_id // M
        metadata["mw"] = mw_cycle[band]
    elif run_id < 8 * M:
        metadata["grid"] = "fine"
        metadata["scenario"] = run_id % M
        band = (run_id - 4 * M) // M
        metadata["mw"] = mw_cycle[band]
    elif run_id == total_non_b0:
        metadata["grid"] = "coarse"
        metadata["mw"] = 0.0
    elif run_id == total_non_b0 + 1:
        metadata["grid"] = "fine"
        metadata["mw"] = 0.0
    else:
        raise ValueError(f"run_id {run_id} exceeds expected range (0-{total_non_b0 + 1})")

    if metadata["mw"] == 0.0:
        metadata["label"] = "B0"
        metadata["rundir"] = f"../geoclaw_output/{metadata['grid']}_B0"
    else:
        metadata["label"] = f"Mw {metadata['mw']:.1f}"
        metadata["rundir"] = (
            f"../geoclaw_output/{metadata['grid']}_{metadata['mw']:.1f}/"
            f"run_{metadata['scenario']}"
        )

    return metadata


def choose_balanced_runs(coarse_ids, fine_ids, n_total):
    """Return a balanced random selection of coarse and fine run IDs."""
    if n_total <= 0:
        return []

    coarse_ids = list(coarse_ids)
    fine_ids = list(fine_ids)

    if not coarse_ids and not fine_ids:
        return []

    n_coarse_target = (n_total + 1) // 2  # coarse receives the extra when odd
    n_fine_target = n_total - n_coarse_target

    selected = []

    if coarse_ids:
        selected.extend(random.sample(coarse_ids, min(n_coarse_target, len(coarse_ids))))
    if fine_ids:
        selected.extend(random.sample(fine_ids, min(n_fine_target, len(fine_ids))))

    remaining_needed = n_total - len(selected)
    if remaining_needed > 0:
        remaining_pool = [
            rid for rid in coarse_ids + fine_ids
            if rid not in selected
        ]
        if remaining_pool:
            selected.extend(random.sample(remaining_pool, min(remaining_needed, len(remaining_pool))))

    return sorted(selected[:n_total])


def chunk_run_ids(run_ids, max_chunks):
    """Split run_ids into at most max_chunks contiguous groups."""
    if max_chunks <= 0:
        return []

    run_ids_sorted = sorted(run_ids)
    if not run_ids_sorted:
        return []

    if max_chunks >= len(run_ids_sorted):
        return [[rid] for rid in run_ids_sorted]

    # Use numpy to split while preserving order and contiguity
    splits = np.array_split(run_ids_sorted, max_chunks)
    return [list(map(int, chunk)) for chunk in splits if len(chunk) > 0]


def parse_cli_args():
    parser = argparse.ArgumentParser(
        description="Run GeoClaw CC CSZ South scenarios in test or full-parallel mode."
    )
    parser.add_argument(
        "--mode",
        choices=["test", "all", "b0"],
        default="test",
        help="Run random test subset ('test'), all scenarios ('all'), or only bathymetry B0 cases ('b0').",
    )
    parser.add_argument(
        "--n-test",
        type=int,
        default=5,
        help="Number of random cases to run in test mode (non-B0 runs only).",
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
        "--include-b0",
        action="store_true",
        help="Include coarse and fine B0 runs in addition to the selected cases.",
    )
    return parser.parse_args()


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
    Total 802 runs:

          run_id + 1
          -----------------------------------------------------------------
          1    100    200    300    400    500    600    700    800 801 802
 grid-size  |           coarse          |           fine            | c | f |
     Mw     |  8.6 |  8.8 |  9.0 |  9.2 |  8.6 |  8.8 |  9.0 |  9.2 | 0.0 0.0

    """

    # helper functions to set slip distribution
    from clawpack.geoclaw import fgmax_tools



    run_id = self._run_id
    etopo_dir = os.path.join(driver_home,"DataFiles")
    topodir = os.path.join(driver_home,"DataFiles")

    # load input info
    if self._input_info == None:
        scn_fname = os.path.join(self._run_home,'scenario_pts.txt') 
        scn = np.loadtxt(scn_fname)
        scn_list = scn.tolist()
    else:
        scn_list = self._input_info
    
    # total number of runs
    M = len(scn_list)
    N = 8*M + 2     # 8*M runs plus two empty bathymetry runs

    if run_id == N:
        raise StopIteration()

    else:
        
        #=========================
        # set coarse and fine grids
        #
        t_shelf = 0.   # time approaching continental slope
        t_harbor = 0.  # time approaching harbor

        if ((run_id >= 0) and (run_id < 4*M)) or (run_id == 8*M):
            #------------------
            # setrun for coarse
            #
            grid = 'coarse'
            self._rundata = setrun_coarse(setgeo_coarse)   # includes fgmax setup

            # empty runs to obtain bathymetry (may be overwritten below)
            dir_grid_Mw = '../geoclaw_output/' + str(grid) + '_B0'
            self._rundir = dir_grid_Mw

            self._rundata.amrdata.amr_levels_max = 4
            # coarse grid run = 10sec
            # dx = 30min, 5min, 1min, 10sec
            self._rundata.amrdata.refinement_ratios_x = [6, 5, 6]
            self._rundata.amrdata.refinement_ratios_y = [6, 5, 6]
            self._rundata.amrdata.refinement_ratios_t = [6, 5, 6]

            # add topography (coarse)
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

        elif ((run_id >= 4*M) and (run_id < 8*M)) or (run_id == 8*M+1):
            #----------------
            # setrun for fine
            #
            grid = 'fine'
            
            self._rundata.amrdata.amr_levels_max = 6

            # fine grid run = 2/3 seconds
            # dx = 30 minutes, 5 minutes, 1 minute, 10 seconds, 2 seconds, 2/3 seconds
            self._rundata.amrdata.refinement_ratios_x = [6, 5, 6, 5, 3]
            self._rundata.amrdata.refinement_ratios_y = [6, 5, 6, 5, 3]
            self._rundata.amrdata.refinement_ratios_t = [6, 5, 6, 5, 3]

            regions = self._rundata.regiondata.regions 
            regions = []  # clear existing entries but keep reference
            # between shelf and CC
            regions.append(\
                    [2, 4, t_shelf, 1e9, -125, -124.05, 40.5, 43]) 
            regions.append(\
                    [4, 5, t_harbor, 1e9, -124.26, -124.14, 41.67,   41.79])
            regions.append(\
                    [6, 6, t_harbor, 1e9, -124.218,-124.17, 41.7345, 41.77])
            self._rundata.regiondata.regions = regions   # ✅ assign back

            # add topography (fine)
            topofiles = self._rundata.topo_data.topofiles
            # for topography, append lines of the form
            #    [topotype, minlevel, maxlevel, t1, t2, fname]
            topofiles[:] = []  # clear existing entries but keep reference

            topofiles.append([3, 1, 6, 0., 1.e10, \
                    os.path.join(etopo_dir, 'etopo1_-130_-124_38_45_1min.asc')])
            topofiles.append([3, 4, 6, 0., 1.e10, \
                    os.path.join(topodir, 'cc-1sec.asc')])
            # topofiles.append([3, 6, 6, 0., 1.e10, \
            #         os.path.join(topodir,'cc-1_3sec-c_pierless.asc')])

        
        # check that topo files exist
        for _, _, _, _, _, f in topofiles:
            if not os.path.exists(f):
                print(f"[FATAL] Topography file missing: {f}")
                raise FileNotFoundError(f)

            
        #
        # set desired magnitude
        #
        if ((run_id >= 0) and (run_id < M)) \
                            or ((run_id >= 4*M) and (run_id < 5*M)):
            self.KL_Mw_desired = 8.6
        elif ((run_id >= M) and (run_id < 2*M)) \
                            or ((run_id >= 5*M) and (run_id < 6*M)):
            self.KL_Mw_desired = 8.8
        elif ((run_id >= 2*M) and (run_id < 3*M)) \
                            or ((run_id >= 6*M) and (run_id < 7*M)):
            self.KL_Mw_desired = 9.0
        elif ((run_id >= 3*M) and (run_id < 4*M)) \
                            or ((run_id >= 7*M) and (run_id < 8*M)):
            self.KL_Mw_desired = 9.2
        
        #
        # set slip distribution
        #
        # run_id_mod = run_id - 100*(run_id/100)
        run_id_mod = run_id % M  # stay within available scenarios
        m = scn_list[run_id_mod]
        self.set_KL_slip(m)
    
        if run_id < 8*M:
            dir_grid_Mw = '../geoclaw_output/' + str(grid) + '_' + str(self.KL_Mw_desired)
            self._rundir = os.path.join(dir_grid_Mw, 'run_' + str(run_id_mod))
        else:
            # empty runs to obtain bathymetry
            
            dir_grid_Mw = '../geoclaw_output/' + str(grid) + '_B0'
            self._rundir = dir_grid_Mw
            self.KL_Mw_desired = 0.0
            self.set_KL_slip([0.]*len(m))   # set output
            self._rundata.clawdata.output_times = [1.0, 3.0]
        
        # --- Compact progress info (after Mw assigned) ---
        if getattr(self, "_progress_enabled", True):
            print(f"[INFO] Case {run_id + 1} / {N} | run_id = {run_id} | Mw = {self.KL_Mw_desired}")
            
        self._run_id += 1
        
        return self


def create_configured_drom():
    """Create a Drom instance with GeoClaw input configured for this scenario."""
    drom = rcrom.Drom()
    drom.GeoClawInput.fault = copy.deepcopy(get_fault())
    drom.GeoClawInput.set_iter(iter_fun)
    drom.GeoClawInput.set_rundata(setrun=setrun_coarse, setgeo=setgeo_coarse)
    drom.GeoClawInput.KL_expand(Lstrike=LSTRIKE, Ldip=LDIP,
                distribution='Lognormal', tau=depth_taper,
                nterms=20, KL_Mw_desired=9.0)
    return drom


def run_case_worker(args):
    """Worker entry-point: run one or more GeoClaw cases identified by run_ids."""
    run_ids, driver_path = args
    run_ids = sorted(run_ids)
    if not run_ids:
        return []

    pending = set(run_ids)
    results = []

    try:
        os.chdir(driver_path)
    except OSError:
        pass

    drom = create_configured_drom()
    drom.GeoClawInput._progress_enabled = False

    for geoinput in drom.GeoClawInput:
        current_id = geoinput._run_id - 1

        if current_id not in pending:
            continue

        run_dir = geoinput._rundir
        Mw = geoinput.KL_Mw_desired

        orig_geo_id = drom.GeoClawInput._run_id
        orig_drom_id = drom._run_id
        drom.GeoClawInput._run_id = current_id
        drom._run_id = current_id

        start = time.perf_counter()
        try:
            print(f"[RUN-{current_id:03d}] Starting GeoClaw in {run_dir} (Mw={Mw})")
            drom.evaluate_hdm()
            duration = time.perf_counter() - start
            print(f"[RUN-{current_id:03d}] Completed successfully in {run_dir} "
                  f"({format_duration(duration)})")
            results.append((current_id, "OK", duration))
        except Exception as exc:
            duration = time.perf_counter() - start
            print(f"[RUN-{current_id:03d}] FAILED in {run_dir}: {exc} "
                  f"({format_duration(duration)})")
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


if __name__=='__main__':
    
    args = parse_cli_args()

    if args.seed is not None:
        random.seed(args.seed)

    # ==========================================================================
    # Build base configuration once (for informational output / caching)
    # ==========================================================================
    get_fault()  # prints fault subdivision information

    scn_path = os.path.join(driver_home, 'scenario_pts.txt')
    try:
        scn = np.loadtxt(scn_path)
    except OSError as err:
        print(f"[ERROR] Failed to read scenario file {scn_path}: {err}")
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
    total_non_b0 = 8 * M
    total_runs = total_non_b0 + 2

    if args.mode == "b0":
        selected_run_ids = []
        if total_runs >= total_non_b0 + 1:
            selected_run_ids.append(total_non_b0)
        if total_runs >= total_non_b0 + 2:
            selected_run_ids.append(total_non_b0 + 1)
        print(f"[INFO] B0 mode: scheduling bathymetry-only runs {selected_run_ids}.")

    elif args.mode == "test":
        if total_non_b0 == 0:
            print("[ERROR] No non-bathymetry runs are available for test mode.")
            sys.exit(1)

        ntest = max(1, min(args.n_test, total_non_b0))
        coarse_ids = list(range(0, 4 * M))
        fine_ids = list(range(4 * M, 8 * M))
        selected_run_ids = choose_balanced_runs(coarse_ids, fine_ids, ntest)

        if len(selected_run_ids) < ntest:
            print(f"[WARN] Only {len(selected_run_ids)} runs available for selection; requested {ntest}.")

        coarse_count = sum(1 for rid in selected_run_ids if rid < 4 * M)
        fine_count = len(selected_run_ids) - coarse_count
        print(f"[INFO] Test mode balanced selection: {coarse_count} coarse, {fine_count} fine runs.")

    else:
        selected_run_ids = list(range(total_runs))
        print(f"[INFO] All mode: scheduling all {total_runs} runs (including coarse/fine B0).")

    if args.include_b0:
        b0_ids = []
        if total_runs >= total_non_b0 + 1:
            b0_ids.append(total_non_b0)
        if total_runs >= total_non_b0 + 2:
            b0_ids.append(total_non_b0 + 1)
        if b0_ids:
            selected_run_ids = sorted(set(selected_run_ids + b0_ids))
            print(f"[INFO] Added B0 runs: {b0_ids}")

    if not selected_run_ids:
        print("[ERROR] No runs selected.")
        sys.exit(1)

    print("[PLAN] Scheduled runs:")
    for rid in selected_run_ids:
        meta = describe_run(rid, M)
        scenario_str = "-" if meta["scenario"] is None else str(meta["scenario"])
        print(f"  - run_id {rid:03d} | grid={meta['grid']:<6} | label={meta['label']:<6} "
              f"| scenario={scenario_str:<3} | rundir={meta['rundir']}")

    cpu_count = mp.cpu_count()
    max_workers = args.processes if args.processes is not None else min(8, cpu_count)
    if max_workers <= 0:
        print("[ERROR] Number of processes must be positive.")
        sys.exit(1)

    run_chunks = chunk_run_ids(selected_run_ids, max_workers)
    if not run_chunks:
        print("[ERROR] Failed to partition run IDs for parallel execution.")
        sys.exit(1)

    actual_workers = len(run_chunks)
    print(f"[INFO] Using {actual_workers} parallel worker(s) (CPU count: {cpu_count}).")

    worker_args = [(chunk, driver_home) for chunk in run_chunks]

    wall_start = time.perf_counter()
    with mp.Pool(processes=actual_workers) as pool:
        worker_results = pool.map(run_case_worker, worker_args)
    wall_total = time.perf_counter() - wall_start

    flat_results = [item for sublist in worker_results for item in sublist]

    ok = sum(1 for _, status, _ in flat_results if status == "OK")
    failures = [(run_id, status, duration) for run_id, status, duration in flat_results if status.startswith("FAIL")]
    skipped = [(run_id, status, duration) for run_id, status, duration in flat_results if status == "SKIPPED"]
    total_run_time = sum(duration for _, _, duration in flat_results)
    total_success_time = sum(duration for _, status, duration in flat_results if status == "OK")

    print("\n[SUMMARY]")
    print(f"  ✅ Successful runs: {ok}")
    if failures:
        for run_id, status, duration in sorted(failures):
            print(f"  ❌ Run {run_id}: {status} ({format_duration(duration)})")
    else:
        print("  ❌ Failed runs: 0")
    if skipped:
        print(f"  [WARN] Skipped runs: {sorted(run_id for run_id, _, _ in skipped)}")
    print(f"  ⏱ Total wall time: {format_duration(wall_total)}")
    print(f"  ⏱ Sum of run durations: {format_duration(total_run_time)}")
    if ok:
        avg_duration = total_success_time / ok
        print(f"  ⏱ Mean successful run: {format_duration(avg_duration)}")
    print("--------------------------------------------------")
