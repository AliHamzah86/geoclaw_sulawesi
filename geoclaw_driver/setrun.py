
import os
import numpy as np

try:
    CLAW = os.environ['CLAW']
except:
    raise Exception("*** Must first set CLAW enviornment variable")

scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')
driver_home = os.getcwd()      # directory where all runs will be done

# Shared data locations (mirror run_scenario.py defaults)
DATA_DIR = os.path.normpath(os.path.join(driver_home, '..', 'DataFiles'))

#------------------------------
def setrun(setgeo,claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data 
    
    
    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"



    t_shelf = 0.   # time approaching continental slope
    t_harbor = 0.  # time approaching harbor

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------


    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim
    
    # Lower and upper edge of computational domain:
    clawdata.lower[0] = -127.5          # xlower
    clawdata.upper[0] = -123.5         # xupper
    clawdata.lower[1] = 38.5          # ylower
    clawdata.upper[1] = 44.5          # yupper
    
    # Number of grid cells:
    clawdata.num_cells[0] =  8      # mx
    clawdata.num_cells[1] = 12      # my

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0
    

    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in 
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00005'  # File to use for restart data
    
    
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
 
    clawdata.output_style = 2
 
    if clawdata.output_style==1:
        # Output ntimes frames at equally spaced times up to tfinal:
        # Can specify num_output_times = 0 for no output
        clawdata.num_output_times = 24
        clawdata.tfinal = 2*3600.
        clawdata.output_t0 = False  # output at initial (or restart) time?
        
    elif clawdata.output_style == 2:
        # Specify a list or numpy array of output times:
        # Include t0 if you want output at the initial time.
        
        # default time and magnitude
        TMID = 600.
        TFINAL = 3600.*2.5
        clawdata.output_times = [float(t) for t in ([2.0] + [float(tt) for tt in np.linspace(TMID, TFINAL, 15)])]


    elif clawdata.output_style == 3:
        # Output every step_interval timesteps over total_steps timesteps:
        clawdata.output_step_interval = 1
        clawdata.total_steps = 1
        clawdata.output_t0 = True  # output at initial (or restart) time?
        

    clawdata.output_format = 'binary'      # 'ascii', 'binary', 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0
    

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 0
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==True:  variable time steps used based on cfl_desired,
    # if dt_variable==Falseixed time steps dt = dt_initial always used.
    clawdata.dt_variable = True
    
    # Initial time step for variable dt.  
    # (If dt_variable==0 then dt=dt_initial for all steps)
    clawdata.dt_initial = 1
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99
    
    # Desired Courant number if variable dt used 
    clawdata.cfl_desired = 0.75
    # max Courant number to allow without retaking step with a smaller dt:
    clawdata.cfl_max = 1.0
    
    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000


    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2
    
    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'
    
    # For unsplit method, transverse_waves can be 
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2
    
    
    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3
    
    # List of limiters to use for each wave family:  
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'vanleer'  ==> van Leer
    #   4 or 'mc'       ==> MC limiter
    clawdata.limiter = ['vanleer', 'vanleer', 'vanleer']
    
    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms
    
    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 1
    
    
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 or 'user'     => user specified (must modify bcNamr.f to use this option)
    #   1 or 'extrap'   => extrapolation (non-reflecting outflow)
    #   2 or 'periodic' => periodic (must specify this at both boundaries)
    #   3 or 'wall'     => solid wall for systems where q(2) is normal velocity
    
    clawdata.bc_lower[0] = 'extrap'   # at xlower
    clawdata.bc_upper[0] = 'extrap'   # at xupper

    clawdata.bc_lower[1] = 'extrap'   # at ylower
    clawdata.bc_upper[1] = 'extrap'   # at yupper
                  
       
    # ---------------
    # Gauges:
    # ---------------

    gauges = rundata.gaugedata.gauges
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]

    # gauges on transects:
    ng = 6
    xg = np.linspace(-124.186, -124.212, ng)
    yg = np.linspace(41.735, 41.76625, ng)
    for j in range(0,ng):
        gaugeno = 201 + j
        gauges.append([gaugeno, xg[j], yg[j], t_harbor, 1.e10])

    ng = 7
    xg = np.linspace(-124.179, -124.205, ng)
    yg = np.linspace(41.735, 41.76625, ng)
    for j in range(0,ng):
        gaugeno = 301 + j
        gauges.append([gaugeno, xg[j], yg[j], t_harbor, 1.e10])


    # horizontal transects to open ocean
    gauges.append([501, -124.3, 41.5, t_shelf, 1.e10])
    gauges.append([502, -124.5, 41.5, t_shelf, 1.e10])
    gauges.append([503, -124.7, 41.5, t_shelf, 1.e10])
    gauges.append([504, -124.9, 41.5, t_shelf, 1.e10])
    gauges.append([505, -125.1, 41.5, t_shelf, 1.e10])
    gauges.append([506, -125.3, 41.5, t_shelf, 1.e10])

    gauges.append([511, -124.3, 41.7, t_shelf, 1.e10])
    gauges.append([512, -124.5, 41.7, t_shelf, 1.e10])
    gauges.append([513, -124.7, 41.7, t_shelf, 1.e10])
    gauges.append([514, -124.9, 41.7, t_shelf, 1.e10])
    gauges.append([515, -125.1, 41.7, t_shelf, 1.e10])
    gauges.append([516, -125.3, 41.7, t_shelf, 1.e10])

    gauges.append([521, -124.3, 42., t_shelf, 1.e10])
    gauges.append([522, -124.5, 42., t_shelf, 1.e10])
    gauges.append([523, -124.7, 42., t_shelf, 1.e10])
    gauges.append([524, -124.9, 42., t_shelf, 1.e10])
    gauges.append([525, -125.1, 42., t_shelf, 1.e10])
    gauges.append([526, -125.3, 42., t_shelf, 1.e10])

    # gauge locations based on relative error
    gauges.append([601, -124.20930028, 41.75447275, t_shelf, 1.e10])
    gauges.append([602, -124.21224899, 41.76223506, t_shelf, 1.e10])
    gauges.append([603, -124.19504245, 41.75001983, t_shelf, 1.e10])
    gauges.append([604, -124.21280513, 41.75523626, t_shelf, 1.e10])
    gauges.append([605, -124.20152458, 41.75949634, t_shelf, 1.e10])
    gauges.append([606, -124.20623252, 41.75453827, t_shelf, 1.e10])
    gauges.append([607, -124.20422907, 41.76285853, t_shelf, 1.e10])
    gauges.append([608, -124.20423643, 41.74837307, t_shelf, 1.e10])
    gauges.append([609, -124.21174351, 41.75235224, t_shelf, 1.e10])
    gauges.append([610, -124.21309088, 41.75877177, t_shelf, 1.e10])

    gauges.append([611, -124.20761466, 41.75697547, t_shelf, 1.e10])
    gauges.append([612, -124.20435671, 41.75754751, t_shelf, 1.e10])
    gauges.append([613, -124.21103431, 41.76499358, t_shelf, 1.e10])
    gauges.append([614, -124.19812103, 41.74777778, t_shelf, 1.e10])
    gauges.append([615, -124.21060714, 41.7572685 , t_shelf, 1.e10])
                        
                  
    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
      # Do not checkpoint at all
      pass

    elif clawdata.checkpt_style == 1:
      # Checkpoint only at tfinal.
      pass

    elif clawdata.checkpt_style == 2:
      # Specify a list of checkpoint times.  
      clawdata.checkpt_times = [0.1,0.15]

    elif clawdata.checkpt_style == 3:
      # Checkpoint every checkpt_interval timesteps (on Level 1)
      # and at the final time.
      clawdata.checkpt_interval = 5

    
    # ---------------
    # AMR parameters:   (written to amr.data)
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = 4

    # List of refinement ratios at each level (length at least amr_level_max-1)

    # grid run = 10"
    # dx = 30', 5', 1', 10"
    amrdata.refinement_ratios_x = [6, 5, 6]
    amrdata.refinement_ratios_y = [6, 5, 6]
    amrdata.refinement_ratios_t = [6, 5, 6]


    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length num_aux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).
    amrdata.aux_type = ['center', 'capacity', 'yleft']


    # Flag for refinement based on Richardson error estimater:
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag_richardson_tol = 1.0  # Richardson tolerance
    
    # Flag for refinement using routine flag2refine:
    amrdata.flag2refine = True      # use this?
    amrdata.flag2refine_tol = 0.5  # tolerance used in this routine
    # Note: in geoclaw the refinement tolerance is set as wave_tolerance below 
    # and flag2refine_tol is unused!

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 10

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width  = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.7

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 0      


    # ---------------
    # Regions:
    # ---------------
    regions = rundata.regiondata.regions 
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    regions.append([1, 2, 0., 1e9, -180, 180, -90, 90])   #whole world
    regions.append([3, 3, 0., 1800, -126, -123.5, 39, 50.5]) #earthquake source 
    regions.append([2, 3, t_shelf, 1e9, -125, -124.05, 40.5, 43]) # between shelf and CC
    regions.append([3, 4, t_harbor, 1e9, -124.26, -124.14, 41.67,   41.79])
    regions.append([4, 4, t_harbor, 1e9, -124.218,-124.17, 41.7345, 41.77])

    #  ----- For developers ----- 
    # Toggle debugging print statements:
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = False      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting
    
    return rundata

    # end of function setrun
# ----------------------

#-------------------
def setgeo(rundata):
#-------------------
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
    refinement_data.speed_tolerance = 0.0
    refinement_data.deep_depth = 100.0
    refinement_data.max_level_deep = 4

    # == settopo.data values ==

    # where to find etopo1 topography:
    etopo_dir = DATA_DIR
    topodir = DATA_DIR   # for other topo files

    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]

    topofiles.append([3, 1, 4, 0., 1.e10, \
            os.path.join(etopo_dir, 'etopo1_-130_-124_38_45_1min.asc')])

    # topofiles.append([3, 3, 4, 0., 1.e10, \
    #         os.path.join(topodir, 'cc-1sec.asc')])

    # == setdtopo.data values ==
    # for moving topography, append lines of the form: [topotype, minlevel,maxlevel,fname]
    rundata.dtopo_data.dtopofiles = []
    dtopofiles = rundata.dtopo_data.dtopofiles
    dtopodir = driver_home
    dtopotype  = 3
    dtopofiles.append([dtopotype, 3, 3, \
            os.path.join(dtopodir,'dtopo.tt3')])

    # == setqinit.data values ==
    rundata.qinit_data.qinit_type =  0
    rundata.qinit_data.qinitfiles = []
    qinitfiles = rundata.qinit_data.qinitfiles 

    from clawpack.geoclaw import fgmax_tools

    fgmax1 = fgmax_tools.FGmaxGrid()
    fgmax1.read_fgmax_grids_data(1, os.path.join(driver_home, 'fgmax1_coarse.txt'))
    fgmax2 = fgmax_tools.FGmaxGrid()
    fgmax2.read_fgmax_grids_data(2, os.path.join(driver_home, 'fgmax2_coarse.txt'))
    fgmax3 = fgmax_tools.FGmaxGrid()
    fgmax3.read_fgmax_grids_data(3, os.path.join(driver_home, 'fgmax3_coarse.txt'))

    fgmax_grids = rundata.fgmax_data.fgmax_grids
    fgmax_grids.append(fgmax1)
    fgmax_grids.append(fgmax2)
    fgmax_grids.append(fgmax3)

    rundata.fgmax_data.num_fgmax_val = 2


    return rundata
    # end of function setgeo
    # ----------------------
