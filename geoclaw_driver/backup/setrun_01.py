"""
Sulawesi PTHA - setrun.py
New-style configuration for Clawpack >= 5.8
- Defines FGmax grids explicitly (fgmax_grids.data)
- Defines FGout grids explicitly (fgout_grids.data)
"""

import os
import numpy as np
from clawpack.clawutil import data
from clawpack.geoclaw import fgmax_tools, fgout_tools

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')
driver_home = os.getcwd()      # directory where all runs will be done


# ===============================================================
def setrun_coarse(setgeo, claw_pkg='geoclaw'):
# ===============================================================
    """
    Define the parameters used for running Clawpack.
    """

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)
    rundata = setgeo_coarse(rundata)

    # ------------------------------------------------------------------
    # Basic Clawpack parameters:
    # ------------------------------------------------------------------
    clawdata = rundata.clawdata
    clawdata.num_dim = num_dim

    # Domain (adjust for Sulawesi region as needed)
    clawdata.lower[0] = 122.0
    clawdata.upper[0] = 126.0
    clawdata.lower[1] = 0.0
    clawdata.upper[1] = 5.0
    clawdata.num_cells[0] = 60
    clawdata.num_cells[1] = 60

    clawdata.num_eqn = 3
    clawdata.num_aux = 3
    clawdata.capa_index = 2

    # Time control
    clawdata.t0 = 0.0
    clawdata.output_style = 1
    clawdata.num_output_times = 24
    clawdata.tfinal = 3.0 * 3600.0
    clawdata.output_t0 = True
    clawdata.output_format = 'binary'

    # Time stepping
    clawdata.dt_variable = True
    clawdata.dt_initial = 0.2
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    clawdata.order = 2
    clawdata.dimensional_split = 'unsplit'
    clawdata.transverse_waves = 2
    clawdata.num_waves = 3
    clawdata.limiter = ['mc', 'mc', 'mc']
    clawdata.use_fwaves = True
    clawdata.source_split = 'godunov'

    # Boundary conditions
    clawdata.num_ghost = 2
    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'
    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # AMR
    amrdata = rundata.amrdata
    amrdata.amr_levels_max = 4
    amrdata.refinement_ratios_x = [4,4,4]
    amrdata.refinement_ratios_y = [4,4,4]
    amrdata.refinement_ratios_t = [4,4,4]
    amrdata.aux_type = ['center', 'capacity', 'yleft']
    amrdata.regrid_interval = 3
    amrdata.regrid_buffer_width = 2
    amrdata.clustering_cutoff = 0.7

    # Regions (example)
    rundata.regiondata.regions = [
        [1, 4, 0., 1e10, 122., 126., 0., 5.]
    ]

    # Gauges (example)
    rundata.gaugedata.gauges = [
        [101, 124.5, 1.0, 0., 1.e10]
    ]

    # ------------------------------------------------------------------
    # GeoClaw-specific parameters
    # ------------------------------------------------------------------
    geo_data = rundata.geo_data
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2
    geo_data.earth_radius = 6367.5e3
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 1e6

    # Topography
    topo_data = rundata.topo_data
    topo_data.topofiles.append([3, 'topo/sulawesi.asc'])

    # Dtopo (example)
    dtopo_data = rundata.dtopo_data
    dtopo_data.dtopofiles.append([3, 'dtopo/dtopo_NS.tt3'])

    # Qinit
    rundata.qinit_data.qinit_type = 0

    # ==========================================================
    # FGmax grids (new style, written to fgmax_grids.data)
    # ==========================================================
    rundata.fgmax_data.num_fgmax_val = 2
    fgmax_grids = rundata.fgmax_data.fgmax_grids

    fg1 = fgmax_tools.FGmaxGrid()
    fg1.fgno = 1
    fg1.point_style = 2
    fg1.x1, fg1.x2 = 122.0, 126.0
    fg1.y1, fg1.y2 = 0.0, 5.0
    fg1.nx, fg1.ny = 50, 50
    fg1.tstart_max = 0.0
    fg1.tend_max = 1.e10
    fg1.dt_check = 60.0
    fg1.min_level_check = 2
    fg1.arrival_tol = 1.e-2
    fgmax_grids.append(fg1)

    # (Add more FGmax grids here if needed)

    # ==========================================================
    # FGout grids (new style, written to fgout_grids.data)
    # ==========================================================
    # fgout_grids = rundata.fgout_data.fgout_grids

    # fgout1 = fgout_tools.FGoutGrid(fgno=1)
    # fgout1.tstart = 0.0
    # fgout1.tend = clawdata.tfinal
    # fgout1.nout = 12
    # fgout1.x1, fgout1.x2 = 122.0, 126.0
    # fgout1.y1, fgout1.y2 = 0.0, 5.0
    # fgout1.nx, fgout1.ny = 60, 60
    # fgout1.output_format = 'binary32'
    # fgout1.output_q_components = [1,2,3]
    # fgout1.output_aux_components = []
    # fgout1.output_aux_onlyonce = True
    # fgout_grids.append(fgout1)

    return rundata

# ===============================================================
def setgeo_coarse(rundata):
# ===============================================================
    """
    Set GeoClaw specific runtime parameters.
    """

    try:
        geo_data = rundata.geo_data
    except:
        print ("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system =  2
    geo_data.earth_radius = 6367500.0

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.     # no tides
    # ******Set in run_tests.py ******

    geo_data.dry_tolerance = 0.001
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 100.0

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.01
    refinement_data.deep_depth = 100.0
    refinement_data.max_level_deep = 4

    # == settopo.data values ==

    # where to find etopo1 topography:
    etopo_dir = driver_home
    topodir = driver_home   # for other topo files
    topodir = os.path.join(topodir, 'DataFiles')

    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]

    topofiles.append([3, 1, 4, 0., 1.e10, \
            os.path.join(etopo_dir, 'etopo1_-130_-124_38_45_1min.asc')])
    topofiles.append([-3, 3, 4, 0., 1.e10, \
            os.path.join(topodir, 'cc-1sec.asc')])

    # == setdtopo.data values ==
    rundata.dtopo_data.dtopofiles = []
    dtopofiles = rundata.dtopo_data.dtopofiles
    # for moving topography, append lines of the form :  
    #   [topotype, minlevel,maxlevel,fname]

    dtopodir = driver_home
    dtopotype  = 3
    dtopofiles.append([dtopotype, 3, 3, \
            os.path.join(dtopodir,'dtopo.tt3')])


    # == setqinit.data values ==
    rundata.qinit_data.qinit_type =  0
    rundata.qinit_data.qinitfiles = []
    qinitfiles = rundata.qinit_data.qinitfiles 
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    # == fixedgrids.data values ==   
    # rundata.fixed_grid_data.fixedgrids = []
    # fixedgrids = rundata.fixed_grid_data.fixedgrids
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    # == fgmax.data values ==
    # fgmax_files = rundata.fgmax_data.fgmax_files
    # fgmax_grids = rundata.fgmax_data.fgmax_grids
    # for fixed grids append to this list names of any fgmax input files
    # fgmax1 = os.path.join(driver_home,'fgmax1_coarse.txt')
    # fgmax2 = os.path.join(driver_home,'fgmax2_coarse.txt')
    # fgmax3 = os.path.join(driver_home,'fgmax3_coarse.txt')

    # fgmax1 = {'fgmax_input_file': os.path.join(driver_home, 'fgmax1_coarse.txt')}
    # fgmax2 = {'fgmax_input_file': os.path.join(driver_home, 'fgmax2_coarse.txt')}
    # fgmax3 = {'fgmax_input_file': os.path.join(driver_home, 'fgmax3_coarse.txt')}

    # fgmax_files.append(fgmax1_fname)  
    # fgmax_files.append(fgmax2_fname)  
    # fgmax_files.append(fgmax3_fname) 

    # fgmax_grids.append(fgmax1)
    # fgmax_grids.append(fgmax2)
    # fgmax_grids.append(fgmax3) 
    
    # rundata.fgmax_data.num_fgmax_val = 2

# -----------------------

    from clawpack.geoclaw import fgmax_tools

    fgmax_grids = rundata.fgmax_data.fgmax_grids

    fgmax1 = fgmax_tools.FGmaxGrid()
    # fgmax1.read_fgmax_txt(os.path.join(driver_home, 'fgmax1_coarse.txt'))
    fgmax1.read_fgmax_grids_data(1, os.path.join(driver_home, 'fgmax1_coarse.txt'))
    fgmax2 = fgmax_tools.FGmaxGrid()
    # fgmax2.read_fgmax_txt(os.path.join(driver_home, 'fgmax2_coarse.txt'))
    fgmax2.read_fgmax_grids_data(2, os.path.join(driver_home, 'fgmax2_coarse.txt'))
    fgmax3 = fgmax_tools.FGmaxGrid()
    # fgmax3.read_fgmax_txt(os.path.join(driver_home, 'fgmax3_coarse.txt'))
    fgmax3.read_fgmax_grids_data(3, os.path.join(driver_home, 'fgmax3_coarse.txt'))

    fgmax_grids.append(fgmax1)
    fgmax_grids.append(fgmax2)
    fgmax_grids.append(fgmax3)

    rundata.fgmax_data.num_fgmax_val = 2

    # just one fgmax grid for testing fixed grid output (October, 2025)
    from clawpack.geoclaw import fgout_tools

    fgout_grids = rundata.fgout_data.fgout_grids

    fgout1 = fgout_tools.FGoutGrid(fgno=1)
    fgout1.tstart = 0.
    fgout1.tend = 1.e10
    fgout1.nout = 10

    fgout1.x1, fgout1.x2 = -125, -124
    fgout1.y1, fgout1.y2 = 41.5, 42.0
    fgout1.nx, fgout1.ny = 50, 50

    fgout1.output_format = 'binary32'   # or 'ascii' or 'binary64'
    fgout1.output_q_components = [1, 2, 3]  # which components to output
    fgout1.output_aux_components = []      # e.g. [1,2]
    fgout1.output_aux_onlyonce = True      # output aux arrays only at t0   

    fgout_grids.append(fgout1)

    return rundata

    # end of function setgeo_coarse
    # ----------------------

    
def setgeo_fine(rundata):
    """
    GeoClaw parameters for fine runs (adds fgmax + fgout)
    """
    fgmax_grids = []
    fgfiles = ['fgmax1_fine.txt', 'fgmax2_fine.txt', 'fgmax3_fine.txt']

    for i, fgfile in enumerate(fgfiles, start=4):   # start numbering at 4
        fg = fgmax_tools.FGmaxGrid()
        fg.read_fgmax_grids_data(i, os.path.join(os.getcwd(), fgfile))
        fg.fgno = i
        fgmax_grids.append(fg)

    rundata.fgmax_data.fgmax_grids = fgmax_grids
    rundata.fgmax_data.num_fgmax_val = 2

    # FGout grid for fine run, unique fgno
    fgout_grids = rundata.fgout_data.fgout_grids
    fgout1 = fgout_tools.FGoutGrid(fgno=101)  # unique ID separate from coarse
    fgout1.tstart = 0.0
    fgout1.tend = 1.e10
    fgout1.nout = 10
    fgout1.x1, fgout1.x2 = -125, -124
    fgout1.y1, fgout1.y2 = 41.5, 42.0
    fgout1.nx, fgout1.ny = 50, 50
    fgout1.output_format = 'binary32'
    fgout_grids.append(fgout1)

    return rundata

    # end of function setgeo_fine
    # ---------------------- 




# ===============================================================
if __name__ == '__main__':
    import sys
    rundata = setrun_coarse(*sys.argv[1:])
    rundata.write()
    print("Created setrun_coarse.py")

    rundata_coarse = setgeo_coarse(rundata)
    rundata_coarse.write()
    print("Created setgeo_coarse.py")

    rundata_fine = setgeo_fine(rundata)
    rundata_fine.write()
    print("Created setgeo_fine.py")