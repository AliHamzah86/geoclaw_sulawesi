"""
Create files such as fine_xyB0_fg1.npy that should be moved to 
all_runs_npy_files directory
"""

from pylab import *
from clawpack.geoclaw import fgmax_tools
import os

resolution_list = ['coarse', 'fine']

geoclaw_output = os.path.abspath('../geoclaw_output')
os.chdir(geoclaw_output)
all_runs_dir = os.path.abspath('all_runs_npy_files')

for resolution in resolution_list:
    B0dir = f'{resolution}_B0'
    print(f"Using files from {B0dir}")

    for fgno in [1, 2, 3]:
        run_fg = fgmax_tools.FGmaxGrid()
        fgmax_grids_path = f'../geoclaw_driver/fgmax_grids_{resolution}.data'
        if not os.path.exists(fgmax_grids_path):
            raise FileNotFoundError(f"Missing fgmax grids data file: {fgmax_grids_path}")
        run_fg.read_fgmax_grids_data(fgno, data_file=fgmax_grids_path)


        outdir = os.path.join(B0dir, '_output')
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
