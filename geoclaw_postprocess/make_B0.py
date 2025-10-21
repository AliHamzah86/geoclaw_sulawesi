"""
Create files such as fine_xyB0_fg1.npy that should be moved to 
all_runs_npy_files directory
"""

from pylab import *
from clawpack.geoclaw import fgmax_tools
import os


def parse_fgmax_legacy(run_fg, fgno, resolution):
    legacy_path = f'../geoclaw_driver/fgmax{fgno}_{resolution}.txt'
    if not os.path.exists(legacy_path):
        raise FileNotFoundError(f"FGmax definition missing: {legacy_path}")

    tokens = []
    detected_fgno = None
    with open(legacy_path) as handle:
        for raw in handle:
            main, _, comment = raw.partition('#')
            main = main.strip()
            comment = comment.strip().lower()
            if 'fgno' in comment and detected_fgno is None:
                detected_fgno = int(float(main)) if main else fgno
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
            run_fg.X = array([xy[0] for xy in coords])
            run_fg.Y = array([xy[1] for xy in coords])
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


def load_fgmax_grid(run_fg, fgno, resolution):
    grids_file = f'../geoclaw_driver/fgmax_grids_{resolution}.data'
    if os.path.exists(grids_file):
        try:
            run_fg.read_fgmax_grids_data(fgno, data_file=grids_file)
            return
        except ValueError:
            pass
    parse_fgmax_legacy(run_fg, fgno, resolution)


resolution_list = ['coarse', 'fine']

geoclaw_output = os.path.abspath('../geoclaw_output')
os.chdir(geoclaw_output)
all_runs_dir = os.path.abspath('all_runs_npy_files')
os.makedirs(all_runs_dir, exist_ok=True)

for resolution in resolution_list:
    B0dir = f'{resolution}_B0'
    print(f"Using files from {B0dir}")

    for fgno in [1, 2, 3]:
        run_fg = fgmax_tools.FGmaxGrid()
        load_fgmax_grid(run_fg, fgno, resolution)

        outdir = os.path.join(B0dir, '_output')
        if not os.path.isdir(outdir):
            print(f"[WARN] Skipping fgno {fgno}: missing directory {outdir}")
            continue
        run_fg.read_output(fgno, outdir)

        X = run_fg.X
        Y = run_fg.Y
        B0 = run_fg.B

        XX = reshape(X, -1, order='F')
        YY = reshape(Y, -1, order='F')
        BB = reshape(B0, -1, order='F')

        xyB0 = array(vstack((XX, YY, BB)).T)
        fname = os.path.join(all_runs_dir, f'{resolution}_xyB0_fg{fgno}.npy')
        save(fname, xyB0)
        print(f"Created {fname}")
