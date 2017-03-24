###############################################################################
from make_haz import hazard_curves_calcv
###############################################################################

import os,sys
from numpy import integer, size
from numpy import array, zeros, linspace, shape, ndim, nan
from numpy import ma, where
from numpy import load, save, loadtxt

###################################                                               #
# This program makes the hazardcurves_probs.npy files for all the subdirectories  #
# specified in the scenario directory, like transect_fg1, transect_fg2, 2Dgrid_fg3#
# The program is presently set up to make these for                               #
#                                                                                 #
#              scenarionames=['scenario_all_pseudoCM_40']                         #
#                                                                                 #
# The USER INPUT below is modified to make choices for a particular job run.      #
#                                                                                 #
# (Note: I just did one scenario per job run, but the program allows for more).   #
# An output description of the job run is stored in output.hc_make_haz_run if     #
# the program is run as follows when hc_make_haz_run.py has been stored in        #
# $FEMA/ptha_tools:                                                               #
#                                                                                 #
# To run: Go to $FEMA/scenarios_1106/scenario_all_pseudoCM_40 and do:             #
#         python $FEMA/ptha_tools/hc_make_haz_run.py > output.hc_make_haz_run     #
####                                                                              #
####                                                                              #
################ USER INPUT:                                                ####### 
run_namedir={}; doubly=False;   #keep doubly=False
if 0:
    #go to the coarse directory for the .npy file for any run number used
    #use the coarse bathymetry 
    type_name='coarse'          #resolution of the scenario, either 'coarse' or 'fine'
    for j in range(100):
        run_namedir[j]='coarse'
if 0:
    #go to the coarse_mod directory for the .npy file for any run number used
    #use the fine bathymetry 
    type_name='fine'            #coarse_mod uses fine bathymetry, so choose 'fine'
    for j in range(100):
        run_namedir[j]='coarse_mod'
if 0:               
    #go to the fine directory for the .npy file for any run number used
    #use the fine bathymetry
    type_name='fine'
    for j in range(100):
        run_namedir[j]='fine'
if 1:
    #go to the pseudo directory for the .npy file for any run number used
    #use the fine bathymetry 
    type_name='fine'
    for j in range(100):
        #run_namedir[j]='pseudo'
        #run_namedir[j]='pseudoCM'
        run_namedir[j]='pseudoCM_40'
        #run_namedir[j]='pseudoCM_dzCC_PE'
        #run_namedir[j]='pseudoCM_dzCC_PE_40'
        #run_namedir[j]='pseudo_40'
if 0:
    #go to the pseudo_dtopo directory for the .npy file for any run number used
    #use the fine bathymetry 
    type_name='fine'
    for j in range(100):
        #run_namedir[j]='pseudo_dtopo'
        run_namedir[j]='pseudo_dtopo_40'

if 0: #EXPERIMENTING, for dtopo, KL probs
    #go to the pseudo_dtopo directory for the .npy file for any run number used
    #use the fine bathymetry, but then overwrite some of the pseudos to go to the
    #fine directory to get them; thereby, getting a total of 40 fine. 20 were the
    #pseudo_dtopo centroids, and the other 20 were 20 nonoverlapping ones in the
    #highest 23 overall based on probability (3 of the 23 were overlapping) 
    #
    doubly=True         #run_namedir is a dictionary, each entry also a dictionary
    type_name='fine'
    pd='pseudo_dtopo'
    for j in range(100):
        run_namedir[j]={86:pd,88:pd,90:pd,92:pd}
    reset=[[43,86],[19,90],[19,88],[19,86],[25,90],[25,88],[25,86],[26,90],[26,88],[26,86],\
           [19,92],[58,90],[58,86],[58,88],[83,90],[83,88],[83,86],[13,86],[13,88],[13,90]]
    for k in range(len(reset)):
        reset_no=reset[k][0]; reset_mag=reset[k][1];
        run_namedir[reset_no][reset_mag]='fine'

if 0: #EXPERIMENTING, for dtopo, uniform probs
    #go to the pseudo_dtopo directory for the .npy file for any run number used
    #use the fine bathymetry, but then overwrite some of the pseudos to go to the
    #fine directory to get them; thereby, getting a total of 40 fine. 20 were the
    #pseudo_dtopo centroids, and the other 20 were 20 nonoverlapping ones in the
    #highest 21 overall based on probability (1 of the 21 were overlapping) 
    #
    doubly=True         #run_namedir is a dictionary, each entry also a dictionary
    type_name='fine'
    pd='pseudo_dtopo'
    for j in range(100):
        run_namedir[j]={86:pd,88:pd,90:pd,92:pd}
    reset=[[99,88],[24,88],[26,88],[27,88],[28,88],[29,88],[30,88],[31,88],[32,88],[33,88],\
           [34,88],[35,88],[36,88],[37,88],[38,88],[39,88],[40,88],[41,88],[42,88],[44,88]]
    for k in range(len(reset)):
        reset_no=reset[k][0]; reset_mag=reset[k][1];
        run_namedir[reset_no][reset_mag]='fine'
if 0:
    #note:  the fine runs, the coarse_mod, and the pseudo runs have resolution fine
    #       fine bathy should be used with the coarse_mod runs and the pseudo runs
    type_name='fine'
    for j in range(100):
        run_namedir[j]='coarse_mod'
    ## Use fine run for the values in reset below
    reset=[43,19,25,26,58,83,13,20,59,99]
    for nr in reset:
        run_namedir[nr]='fine'
#
#########  choose whether old or new files (new has different magnitudes)         #
#
if 1:                #choose this for FEMA 1106 work
    old_new='1106'
if 0:
    old_new='new'
if 0:
    old_new='old'

#########  Pick one of the following for this job run             #
######### SCENARIOS for the 400 fine and 400 coarse data set      #
#
#scenarionames=['scenario_all_fine']
#scenarionames=['scenario_all_fine_uniform']
#scenarionames=['scenario_all_pseudo']
#scenarionames=['scenario_all_pseudoCM']
scenarionames=['scenario_all_pseudoCM_40']
#scenarionames=['scenario_all_pseudoCM_dzCC_PE']
#scenarionames=['scenario_all_pseudoCM_dzCC_PE_40']
#scenarionames=['scenario_all_pseudo_40']
#scenarionames=['scenario_all_pseudo_dtopo']
#scenarionames=['scenario_all_pseudo_dtopo_40']
#scenarionames=['scenario_all_pseudo_dtopo_addto']
#scenarionames=['scenario_all_pseudo_dtopo_uniform']
#scenarionames=['scenario_all_pseudo_dtopo_uniform_addto']
#scenarionames=['scenario_all_pseudo_uniform']
#scenarionames=['scenario_all_coarse']
#scenarionames=['scenario_all_coarse_mod']
#scenarionames=['scenario_highest10_each']
#scenarionames=['scenario_dzCC_dtopo_etamax_20cls']
#scenarionames=['scenario_dzCC_PE_20cls']
#scenarionames=['scenario_dzCC_PE_40cls']
#######
#scenarionames=['scenario_dtopo_40cls']
#scenarionames=['scenario_coarse_eta_40cls']
#######
#scenarionames=['scenario_coarse_eta_20cls']
#scenarionames=['scenario_coarse_etaCM_20cls']
#scenarionames=['scenario_coarse_etaCM_40cls']
#
#scenarionames=['scenario_highest20_overall']
#scenarionames=['scenario_highest20_uniform']
#scenarionames=['scenario_c_eta_20cls_uniform']
#scenarionames=['scenario_dzCC_dtopo_20cls_uniform']
#
############# All scenarios below here were testing or with OLD 100 data sets ###
#scenarionames=['scenario_all_fine_uniform']
#scenarionames=['scenario_uniform10']
#scenarionames=['scenario_highest10']
#scenarionames=['scenario_highest10_newzetas']
#scenarionames=['scenario_all_fine']
#scenarionames=['scenario_all_fine_newzetas']
#scenarionames=['scenario_all_coarse']
#scenarionames=['scenario_all_coarse_newzetas']
#scenarionames=['scenario_all_coarse_mod']
#scenarionames=['scenario_h10_coarse_mod']

#
#scenario_all_fine is FEMA project, fine runs_0 thru runs_99, 2 transects, 1 2Dgrid
#scenario_all_coarse is FEMA project, coarse runs_0 thru runs_99, 2 transects, 1 2Dgrid
#scenario_highest10  is FEMA project, the 10 fine runs with highest cond. probs
#                                   , 2 transects, 1 2Dgrid
#scenario_uniform10  is FEMA project, the same 10 fine runs above but with
#                                   , uniform conditional probabilities set to .1
#                                   , 2 transects, 1 2Dgrid
#
#scenario_all_fine_uniform is FEMA project, all fine runs with uniform probabilities
#                             , of .01

###  USER SETS zeta exceedance values, using 1261 for the FEMA 1106 results
zeta1=list(linspace(0,12,1201)); zeta2=list(linspace(12.2,24,60))
zeta = array(zeta1+zeta2)
print 'The 1261 exceedance values were zeta=linspace(0,12,1201) + linspace(12.2,24,60)'
print ' '
#############    END OF USER INPUT                    ################################
###########################DONT CHANGE BELOW THIS LINE################################

projectdir = os.environ['FEMA']
print 'INSIDE DRIVER hc_make_haz_run.py, data echo: '
print ' '
print 'The project directory was: ',projectdir
print ' '

for scenarioname in scenarionames:
    scenariodir = projectdir + '/scenarios_' + old_new +'/' + scenarioname
    scenario_weights_file=scenariodir + '/scenario_prb_wgts.txt'
    runs_todo_file = scenariodir + '/runs_todo.txt'
    fg_types_todo_file = scenariodir + '/fg_types_todo.txt'
    fg_transects_todo_file = scenariodir + '/fg_transects_todo.txt'
    fg_2Dgrids_todo_file = scenariodir + '/fg_2Dgrids_todo.txt'

    print 'The scenario directory was: ',scenariodir
    print ' '
    
    runs_todo = loadtxt(runs_todo_file,dtype=integer)
    sh=shape(runs_todo)
    num_runs=sh[0]

    fg_types_todo = loadtxt(fg_types_todo_file,dtype=integer)
    print 'fg_types_todo where 1 means transects, 2 means 2Dgrids is: '
    print fg_types_todo
    print ' '

    events_directories=[]
    if (ndim(runs_todo)==1):   
        for event in runs_todo:
            runsdir=projectdir + '/all_runs_npy_files_old/'
            run_name=run_namedir[event]
            events_directories.append(runsdir + run_name +'_runs' + '/run_'+str(event))
        print 'The master run numbers used for hazard curves are: '
    else:
        for k in range(num_runs): 
            event=runs_todo[k,0]; mag=runs_todo[k,1];
            #####
            if 0:
                runsdir=projectdir + '/all_runs_npy_files_new/'
            if 1:
                runsdir=projectdir + '/redclaw_1106/all_runs_npy_files/'
            #####
            if (doubly==True):
                run_name=run_namedir[event][mag]
            if (doubly==False):
                run_name=run_namedir[event]
            mag=mag/10.0
            strmag='_' + str(mag)
            events_directories.append(runsdir+run_name+strmag+'_runs'+'/run_'+str(event))
        print 'The master run numbers and magnitudes used for hazard curves are: '
    print runs_todo
    print ' '

    print 'The directories for all the events was: '
    print events_directories 
    print ' '

    events_probs = loadtxt(scenariodir + '/scenario_prb_wgts.txt')
    print 'The conditional probability weights for scenario events was: '
    print events_probs
    print ' '

    ######  Now loop over all the work to be done, and call hazard_curves_calcv ###
    ######  for each fixed grid of each type.
    for itype in fg_types_todo:
        if (itype == 1):            #Do all the transect type of fixed grids
            fg_transects_todo=loadtxt(fg_transects_todo_file,dtype=integer)
            if (size(fg_transects_todo) == 1):
                fg_transects_todo = array([fg_transects_todo])
            for iwhich in fg_transects_todo:
                outdir = scenariodir + '/' + 'transect_fg' + str(iwhich)
                hazard_curves_file = outdir + '/hazardcurves_probs.npy'
                fixed_grid_file = runsdir+ type_name +\
                                  '_xyB0_fg'+str(iwhich) +'.npy'

                print '#####################'
                print 'CALLING hazard_curves_calcv with fixed grid number: ',iwhich
                print 'The fixed grid xyB0 file was: ',fixed_grid_file
                print ' '
                print 'The hazard curves .npy file was: '
                print hazard_curves_file
                print ' '
                hazard_curves_calcv(hazard_curves_file,fixed_grid_file,\
                                    events_directories,events_probs,\
                                    zeta,iwhich,comments_file=None)
                print ' '
        if (itype == 2):            #Do all the 2D grid type of fixed grids
            fg_2Dgrids_todo=loadtxt(fg_2Dgrids_todo_file,dtype=integer)
            if (size(fg_2Dgrids_todo) == 1):
                fg_2Dgrids_todo = array([fg_2Dgrids_todo])
            for iwhich in fg_2Dgrids_todo:
                outdir = scenariodir + '/' + '2Dgrid_fg' + str(iwhich)
                hazard_curves_file = outdir + '/hazardcurves_probs.npy'
                fixed_grid_file = runsdir+ type_name +\
                                     '_xyB0_fg'+str(iwhich) +'.npy'

                print '#####################'
                print 'CALLING hazard_curves_calcv with fixed grid number: ',iwhich
                print 'The fixed grid xyB0 file was: ',fixed_grid_file
                print ' '
                print 'The hazard curves .npy file was: '
                print hazard_curves_file
                print ' '
                hazard_curves_calcv(hazard_curves_file,fixed_grid_file,\
                                    events_directories,events_probs,\
                                    zeta,iwhich,comments_file=None)
                print ' '
