from hazard_curves_2Dgrid_plot import hazard_curves_2Dgrid_plot
from hazard_curves_transect_plot import hazard_curves_transect_plot

from numpy import array,save,loadtxt,linspace
from numpy import integer,size

import matplotlib.pyplot as plt

import os

#### The program produces plots for the scenarios described in scenarionames. This#
#### includes plots for all the subdirectories of the scenario, such as           #
#### transect_fg1, transect_fg2, and 2Dgrid_fg3.  I executed it one scenario per  #
#### I executed it one scenario per job run, but a list of scenarios can be done. #
#### It also creates and stores the ZetaFloods.npy files for each fixed grid for  #
#### each scenario.  This file is used when scenarios are compared by the program #
#### hc_compare_run.py.                                                           #                        
####                                                                              #
#### In addition to the plots, an output file is stored in each of the subdirs    #
#### called hc_plot_run.output with specifics about the plots.  Also an output    #
#### file of the job run history called output.hc_plot_run is stored in the       #
#### scenario directory if the program is run as follows.  For example,           #
####                                                                              #
#### To run:  Go to $FEMA/scenarios_1106/scenario_all_SVD_40 and do:              #
####          python $FEMA/ptha_tools/hc_plot_run.py > output.hc_plot_run         #
####                                                                              #
####          when all python programs are stored in $FEMA/ptha_tools             #
####
#### The program calls hazard_curves_2Dgrid_plot and hazard_curves_transect_plot. #
###################################################################################
################ USER INPUT:                                         ##############
#                                                                                 #
if 0:
    type_name='coarse'     #coarse run choose this one
if 1:
    type_name='fine'       #coarse_mod, pseudo, fine all choose this one

if 1:                      #always choose for FEMA 1106 results
    old_new='1106'
if 0:
    old_new='old'
if 0:
    old_new='new'
#
####  Choose which scenarionames we are making plots for.  I did one per job run,
####  but the program can handle a list of them.
####
####  Scenarios for the NEW 400 data coarse and 400 data fine sets #
#scenarionames=['scenario_all_fine']
#scenarionames=['scenario_all_fine_uniform']
#scenarionames=['scenario_all_pseudo']
#scenarionames=['scenario_all_pseudoCM']
#scenarionames=['scenario_all_pseudoCM_40']
#
#DonSub
scenarionames=['scenario_all_SVD_40']
#
#scenarionames=['scenario_all_pseudoCM_dzCC_PE']
#scenarionames=['scenario_all_pseudoCM_dzCC_PE_40']
#scenarionames=['scenario_all_pseudo_40']
#scenarionames=['scenario_all_pseudo_dtopo']
#scenarionames=['scenario_all_pseudo_dtopo_40']
#scenarionames=['scenario_all_pseudo_dtopo_addto']
#scenarionames=['scenario_all_pseudo_dtopo_uniform_addto']
#scenarionames=['scenario_all_pseudo_dtopo_uniform']
#scenarionames=['scenario_all_pseudo_uniform']
#scenarionames=['scenario_all_coarse']
#scenarionames=['scenario_all_coarse_mod']
#scenarionames=['scenario_highest10_each']
#######
#scenarionames=['scenario_dtopo_40cls']
#scenarionames=['scenario_coarse_eta_40cls']
#######
#scenarionames=['scenario_dzCC_dtopo_etamax_20cls']
#scenarionames=['scenario_dzCC_PE_20cls']
#scenarionames=['scenario_dzCC_PE_40cls']
#scenarionames=['scenario_coarse_eta_20cls']
#scenarionames=['scenario_coarse_etaCM_20cls']
#scenarionames=['scenario_coarse_etaCM_40cls']
#scenarionames=['scenario_highest20_overall']
#scenarionames=['scenario_highest20_uniform']
#scenarionames=['scenario_c_eta_20cls_uniform']
#scenarionames=['scenario_dzCC_dtopo_20cls_uniform']
#
#################### Scenarios below are for the OLD 100 data set or tests ########
#scenarionames=['scenario_all_fine']
#scenarionames=['scenario_all_fine_newzetas']
#scenarionames=['scenario_highest10']
#scenarionames=['scenario_highest10_newzetas']
#scenarionames=['scenario_uniform10']
#scenarionames=['scenario_all_fine_uniform']
#scenarionames=['scenario_all_coarse']
#scenarionames=['scenario_all_coarse_newzetas']
#scenarionames=['scenario_all_coarse_mod']  #choose type_name='fine' with this
#scenarionames=['scenario_h10_coarse_mod']  #choose type_name='fine' with this
#                                                                                 #
########END OF USER INPUT     #####################################################

projectdir = os.environ['FEMA']
print 'INSIDE DRIVER hc_plot_run.py, data echo: '
print ' '
print 'The project directory was: ',projectdir
print ' '

#set the directory where the fixed grid files are
if 0:
    runsdir=projectdir + '/all_runs_npy_files_' + old_new + '/'

if 1: #Always use for the FEMA 1106 results
    runsdir=projectdir + '/redclaw_1106/all_runs_npy_files'+ '/'

#Read in the file below only once, since same one used for all transects and 2Dgrids
fixed_big_grid_file = projectdir + '/DataFiles/' + 'CCtopo_1_3sec.tt3'
print 'The fixed grid file for the biggest fixed grid region used for plotting '
print 'little map in upper right graph corner: '
print fixed_big_grid_file
print ' '
from clawpack.geoclaw.topotools import Topography
CCtopo = Topography()
CCtopo.read(fixed_big_grid_file,3)

for scenarioname in scenarionames:
    scenariodir = projectdir + '/scenarios_' + old_new + '/' + scenarioname
    fg_types_todo_file = scenariodir + '/fg_types_todo.txt'
    fg_transects_todo_file = scenariodir + '/fg_transects_todo.txt'
    fg_2Dgrids_todo_file = scenariodir + '/fg_2Dgrids_todo.txt'

    print 'The scenario directory was: ',scenariodir
    print ' '
    
    fg_types_todo = loadtxt(fg_types_todo_file,dtype=integer)
    if (size(fg_types_todo) == 1):
        fg_types_todo = array([fg_types_todo])
    print 'fg_types_todo where 1 means transects, 2 means 2Dgrids is: '
    print fg_types_todo
    print ' '

    depth=True; speed=False; mflux=False;

    if (depth == True):  #using depth only now
        #CHANGE zeta and whichzetas here if desired
        zeta1=list(linspace(0,12,1201)); zeta2=list(linspace(12.2,24,60))
        zeta = array(zeta1+zeta2)
        whichzetas=[0,100,200,300,400,600,800,1000,1200]  #zeta=0.0,1.0,2.0,3.0,4.0,6.0
                                                          #     8.0,10.0,12.0 meters
                                                          #for p-contours or p-transects
        whichzetas=array(whichzetas)

    elif (speed==True): #fix later if used
        zeta1=[0.0,.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5]
        zeta2=[10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0]
        zeta3=[14.5,15.0,15.5,16.0,16.5,17.0]
        zeta=array(zeta1+zeta2+zeta3)
        whichzetas=[0,2,4,6,8,10,12,14,16,18,20]  #zeta=0.0,1,2,3,4,5,6,7,8,9,10 m/sec.
        whichzetas=array(whichzetas)

    else: #mflux is true, fix later if used
        zeta=[0.01, 0.05] + list(linspace(0.1,2,20)) + [4,6,8,10,20,50,100,200,400,600,1000,1500,1700]
        zeta=array(zeta)
        whichzetas=[2,11,25,28]
        whichzetas=array(whichzetas)

    print 'The whichzetas for the p-transects or p-contours were: '
    print whichzetas
    print ' '

    pbarlist=[.1,.2,.3,.4,.5,.6,.7,.8,.9,.95]   #10 zeta-contours or 10 zeta-transects
    print 'The pbarlist for the 10 zeta-contours or 10 zeta-transects was: '
    print pbarlist
    print ' '

    ######  Now loop over all the work to be done for this particular scenario, and ###
    ######  call hazard_curves_transect_plot or hazard_curves_2Dgrid_plot for each  ###
    ######  fixed grid of the appropriate type.                                     ###
    for itype in fg_types_todo:
        if (itype == 1):            #Do all the transect type of fixed grids
            fg_transects_todo=loadtxt(fg_transects_todo_file,dtype=integer)

            if (size(fg_transects_todo) == 1):
                fg_transects_todo = array([fg_transects_todo])

            for iwhich in fg_transects_todo:
                outdir = scenariodir + '/transect_fg' + str(iwhich)
                hazard_curves_file = outdir + '/hazardcurves_probs.npy'
                fixed_grid_file = runsdir+type_name+\
                                  '_xyB0_fg'+str(iwhich) +'.npy'

                #The plotter will automatically choose evenly distributed points for
                #the hazard curves.  If these are good enough, choose whichrows=None;
                #otherwise, make whichrows a list (say of gauge indices) and the plotter
                #will append the evenly distributed points to the whichrows list.
                #BUT NOTE THAT the gauge indices can be transect dependent, and if this
                #is a priority, we would read them in from right scenario directory,
                #and make_scenario.py would have to be modified to put them there.
             
                whichrows=None     

                print '#####################'
                print 'CALLING hazard_curves_transect_plot for fixed grid number: ',iwhich
                print 'The fixed grid xyB0 file was: ',fixed_grid_file
                print ' '
                print 'The hazard curves .npy file was: '
                print hazard_curves_file
                print ' '

                zetapbar=hazard_curves_transect_plot(hazard_curves_file,fixed_grid_file,\
                          zeta,whichzetas,outdir,pbarlist,\
                          CCtopo.X,CCtopo.Y,CCtopo.Z,\
                          whichrows=whichrows,speed=speed,depth=depth,mflux=mflux)

                #Write zetapbar to a file.
                if 1:
                    if (depth==True):
                        outZetafile=outdir + '/' + 'ZetaFloods.npy'
                    elif (speed==True):
                        outZetafile=outdir + '/' + 'SpeedFloods.npy'
                    else:
                        outZetafile=outdir + '/' + 'MfluxFloods.npy'
                    save(outZetafile,zetapbar)
                print ' '

        if (itype == 2):            #Do all the 2D grid type of fixed grids
            fg_2Dgrids_todo=loadtxt(fg_2Dgrids_todo_file,dtype=integer)
            if (size(fg_2Dgrids_todo) == 1):
                fg_2Dgrids_todo = array([fg_2Dgrids_todo])
            for iwhich in fg_2Dgrids_todo:
                outdir = scenariodir +  '/2Dgrid_fg' + str(iwhich)
                hazard_curves_file = outdir + '/hazardcurves_probs.npy'
                fixed_grid_file = runsdir+type_name+\
                                  '_xyB0_fg'+str(iwhich) +'.npy'

                ###      Read Npts1 and Npts2 from file region_dimensions.txt in outdir that
                ###      make_scenerio.py put there. The file will be of length 2 with
                ###      Npts1 and Npts2.
                dimension_file=outdir + '/' + 'region_dimensions.txt'
                dimensions=loadtxt(dimension_file,dtype=integer)
                Npts1=dimensions[0]; Npts2=dimensions[1];

                #The plotter will automatically choose evenly distributed points for
                #the hazard curves.  If these are good enough, choose whichrows=None;
                #otherwise, make whichrows a list (say of gauge indices) and the plotter
                #will append the evenly distributed points to the whichrows list.
                #BUT NOTE THAT the gauge indices will differ between 2Dgrids, and if this
                #is a priority, we would read them in from right scenario directory,
                #and make_scenario.py would have to be modified to put them there.

                whichrows=None

                print '#####################'
                print 'CALLING hazard_curves_2Dgrid_plot with fixed grid number: ',iwhich
                print 'The fixed grid xyB0 file was: ',fixed_grid_file
                print ' '

                print 'The 2Dgrid was Npts1 x Npts2 points in longxlat dimensions'
                print Npts1,Npts2
                print ' '

                print 'The hazard curves .npy file was: '
                print hazard_curves_file
                print ' '
                zetapbar=hazard_curves_2Dgrid_plot(hazard_curves_file,fixed_grid_file,\
                              zeta,whichzetas,outdir,pbarlist,Npts1,Npts2,\
                              CCtopo.X,CCtopo.Y,CCtopo.Z,\
                              whichrows=whichrows,speed=speed,depth=depth,mflux=mflux)
                #Write zetapbar to a file.
                if 1:
                    if (depth==True):
                        outZetafile=outdir + '/' + 'ZetaFloods.npy'
                    elif (speed==True):
                        outZetafile=outdir + '/' + 'SpeedFloods.npy'
                    else:
                        outZetafile=outdir + '/' + 'MfluxFloods.npy'
                    save(outZetafile,zetapbar)
                print ' '
