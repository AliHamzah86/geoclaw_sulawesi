from hazard_curves_2Dgrid_compare import hazard_curves_2Dgrid_compare
from hazard_curves_transect_compare import hazard_curves_transect_compare

from numpy import array,loadtxt,linspace
from numpy import integer,size

import matplotlib.pyplot as plt

import os,sys

###################################################################################
################ USER CHOOSES THE SCENARIOS FOR WHICH COMPARISONS ARE DONE ########
################ AND SETS THE NUMBER OF THE COMPARISON DIRECTORY FOR RESULTS ######
####
####       The program is presently set up to do comparison_number=43        ######
####       The directory comparison_43 and its subdirectory 2Dgrid_fg3       ######
####       were created manually before this program was executed.  Some     ######
####       comparisons also require the subdirectories transect_fg1 and      ######
####       transect_fg2 to also be created before execution of this pgm.     ######
####
####       The program accesses the scenarios directories of all those       ######
####       being compared, and produces comparison plots and comments in     ######
####       the proper comparisondir below. Overall job run information is    ######
####       stored in $FEMA/comparisons_1106/comparison_43 in the file        ######
####       output.hc_compare_run when the program is run as follows:         ######
####
####     To run:  Go to $FEMA/comparisons_1106/comparison_43  and do:        ######
####       python $FEMA/ptha_tools/hc_compare_run.py > output.hc_compare_run ######  
####                                                                         ######
####       when all python programs are stored in $FEMA/ptha_tools           ######                     
####                                                                         ######
####       For completeness, the description of all our comparisons are left ######
####       below for a history of what we did in this project                ######
####                                                                         ######
####     The user will make choice by editing the program below.             ######

if 1:  #Always choose this for the FEMA 1106 report results.
    old_new='1106'
if 0:
    old_new='new'
if 0:
    old_new='old'
#
############## Comparisons for the NEW 400 coarse and 400 fine data set           #
#comparison_number=1
#scenarionames=['scenario_all_fine','scenario_highest10_each']
#scenarioresolution=['fine','fine']

#comparison_number=40
#scenarionames=['scenario_all_fine','scenario_highest10_each',\
#               'scenario_dtopo_40cls','scenario_coarse_eta_40cls']
#scenarioresolution=['fine','fine','fine','fine']

#comparison_number=41
#scenarionames=['scenario_all_fine','scenario_highest10_each',\
#               'scenario_dtopo_40cls','scenario_coarse_eta_40cls',\
#               'scenario_all_pseudo_dtopo_40','scenario_all_pseudo_40']
#scenarioresolution=['fine','fine','fine','fine','fine','fine']

#comparison_number=42
#scenario_dzCC_PE_40cls used CM
#scenarionames=['scenario_all_fine','scenario_highest10_each',\
#               'scenario_dzCC_PE_40cls','scenario_coarse_etaCM_40cls',\
#               'scenario_all_pseudoCM_dzCC_PE_40','scenario_all_pseudoCM_40']
#scenarioresolution=['fine','fine','fine','fine','fine','fine']

comparison_number=43
scenarionames=['scenario_all_fine','scenario_highest10_each',\
               'scenario_all_pseudoCM_40','scenario_all_SVD_40']
scenarioresolution=['fine','fine','fine','fine']

#comparison_number=13
#scenarionames=['scenario_all_fine','scenario_highest10_each',\
#               'scenario_all_pseudo_dtopo','scenario_all_pseudo_dtopo_addto']
#scenarioresolution=['fine','fine','fine','fine']

#comparison_number=14
#scenarionames=['scenario_all_fine_uniform','scenario_dzCC_dtopo_20cls_uniform',\
#               'scenario_all_pseudo_dtopo_uniform','scenario_all_pseudo_dtopo_uniform_addto']
#scenarioresolution=['fine','fine','fine','fine']

#comparison_number=2
#scenarionames=['scenario_all_fine','scenario_all_coarse']
#scenarioresolution=['fine','coarse']

#comparison_number=3
#scenarionames=['scenario_all_fine','scenario_highest10_each','scenario_all_coarse']
#scenarioresolution=['fine','fine','coarse']
#
#comparison_number=4
#scenarionames=['scenario_all_fine','scenario_highest20_overall',\
#               'scenario_dzCC_dtopo_etamax_20cls','scenario_coarse_eta_20cls']
#scenarioresolution=['fine','fine','fine','fine']
#
#comparison_number=5
#scenarionames=['scenario_all_fine','scenario_all_pseudo',\
#               'scenario_dzCC_dtopo_etamax_20cls','scenario_coarse_eta_20cls',\
#               'scenario_all_pseudo_dtopo']
#scenarioresolution=['fine','fine','fine','fine','fine']
#

#comparison_number=15
#scenarionames=['scenario_all_fine','scenario_all_coarse_mod',\
#               'scenario_all_pseudoCM_dzCC_PE','scenario_all_pseudoCM',\
#               'scenario_all_pseudo_dtopo','scenario_all_pseudo']
#scenarioresolution=['fine','fine','fine','fine','fine','fine']

#comparison_number=16
#scenario_dzCC_PE_20cls used CM
#scenarionames=['scenario_all_fine','scenario_all_coarse_mod',\
#                'scenario_dzCC_dtopo_etamax_20cls','scenario_coarse_eta_20cls',\
#                'scenario_dzCC_PE_20cls','scenario_coarse_etaCM_20cls']
#scenarioresolution=['fine','fine','fine','fine','fine','fine']

#comparison_number=6
#scenarionames=['scenario_all_fine_uniform','scenario_all_pseudo_uniform',\
#               'scenario_c_eta_20cls_uniform','scenario_highest20_uniform']
#scenarioresolution=['fine','fine','fine','fine']
#comparison_number=7
#scenarionames=['scenario_all_fine_uniform','scenario_all_pseudo_uniform',\
#               'scenario_dzCC_dtopo_20cls_uniform',\
#               'scenario_c_eta_20cls_uniform',\
#               'scenario_all_pseudo_dtopo_uniform']
#scenarioresolution=['fine','fine','fine','fine','fine']

#comparison_number=10
#scenarionames=['scenario_all_fine','scenario_highest20_overall']
#scenarioresolution=['fine','fine']
#
#comparison_number=11
#scenarionames=['scenario_all_fine_uniform','scenario_highest20_uniform']
#scenarioresolution=['fine','fine']

### Answering the question of whether coarse or coarse_mod alone without any pseudo-fine
### is competitive. Probably concluding both pseudo strategies (based on two different
### clustering techniques is better than the coarse or coarse_mod alone. Not bothering
### with the uniform distribution, since haven't run all_coarse_uniform or all_coarse_mod_uniform
### since can probably rule out with this test alone.
#comparison_number=12
#scenarionames=['scenario_all_fine','scenario_all_coarse',\
#               'scenario_all_pseudo','scenario_all_coarse_mod',\
#               'scenario_all_pseudo_dtopo']
#scenarioresolution=['fine','coarse','fine','fine','fine']
##############

############## comparisons below were for the OLD 100 data set or were tests
#comparison_number=1
#scenarionames=['scenario_all_fine_uniform','scenario_uniform10']
#scenarioresolution=['fine','fine']

#comparison_number=2
#scenarionames=['scenario_all_fine_uniform','scenario_uniform10',\
#              'scenario_all_fine','scenario_highest10']
#scenarioresolution=['fine','fine','fine','fine']

#comparison_number=5
#scenarionames=['scenario_all_fine_newzetas','scenario_all_coarse_mod','scenario_all_coarse_newzetas']
#scenarioresolution=['fine','fine','coarse']

#comparison_number=6
#scenarionames=['scenario_all_fine_newzetas','scenario_highest10_newzetas']
#scenarioresolution=['fine','fine']

#comparison_number=7
#scenarionames=['scenario_all_fine_newzetas','scenario_highest10_newzetas','scenario_h10_coarse_mod',\
#               'scenario_all_coarse_mod']
#scenarioresolution=['fine','fine','fine','fine']


######  END OF USER INPUT                              ##########
############################Dont change below this line##########

projectdir = os.environ['FEMA']
print 'INSIDE DRIVER hc_compare_run.py, data echo: '
print ' '
print 'The project directory was: ',projectdir
print ' '

if 0:
    runsdir=projectdir + '/all_runs_npy_files_' + old_new + '/'
if 1:
    runsdir=projectdir + '/redclaw_1106/all_runs_npy_files' +'/'

#Read in the file below only once, since same one used for all transects and 2Dgrids
fixed_big_grid_file = projectdir + '/DataFiles/' + 'CCtopo_1_3sec.tt3'
print 'The fixed grid file for the biggest fixed grid region used for plotting '
print 'little map in upper right graph corner: '
print fixed_big_grid_file
print ' '
from clawpack.geoclaw.topotools import Topography
CCtopo = Topography()
CCtopo.read(fixed_big_grid_file,3)

scenariodir=[]               #going to be a list
fg_types_todo=[]             #going to be a list of arrays
j=0
for scenarioname in scenarionames:
    scdir=projectdir + '/scenarios_' + old_new + '/' + scenarioname
    scenariodir.append(scdir)
    fg_types_todo_file = scdir + '/fg_types_todo.txt'
    types_array=loadtxt(fg_types_todo_file,dtype=integer)
    if (size(types_array) == 1):
        types_array = array([types_array])
    fg_types_todo.append(types_array)
    j+=1

print 'the scenariodir file was: '
print scenariodir
print ' '
print 'the fg_types_todo list was: '
print fg_types_todo
print ' '

fg_transects_todo_file=[]
fg_2Dgrids_todo_file=[]
num_scenarios = len(scenarionames)
for j in range(num_scenarios):
    if (set(list(fg_types_todo[0])) == set(list(fg_types_todo[j]))): #scenarios same types
        for itype in fg_types_todo[0]:
            if (itype == 1):
                fg_transects_todo_file.append(scenariodir[j] + '/fg_transects_todo.txt')
            if (itype == 2):
                fg_2Dgrids_todo_file.append(scenariodir[j] + '/fg_2Dgrids_todo.txt')
    else:
        print 'can not compare these scenarios. Their todo files have different types'
        print 'scenario number j different type from scenario number 0, j= ',j
        print 'fg_types_todo[0] was: ',fg_types_todo[0]
        print 'fg_types_todo[j] was: ',fg_types_todo[j]
        sys.exit(1)

print 'the fg_transects_todo_file list was: '
print fg_transects_todo_file
print ' '
print 'the fg_2Dgrids_todo_file list was: '
print fg_2Dgrids_todo_file
print ' '

######  Now loop over all the work that has been done for these scenarios, and ###
######  check which transects and 2D grids to compare                          ###
transect_found=0; found=0;
for itype in fg_types_todo[0]:    #set(list(fg_types_todo[j]))=set(list(fg_types_todo[0]))##
    if (itype == 1):              #Do all the transect type of fixed grids
        transect_found=1
        fg_transects_todo=[]      #A list of arrays
        for j in range(num_scenarios):
            fg_transects_todo.append(loadtxt(fg_transects_todo_file[j],dtype=integer))
            if (size(fg_transects_todo[j]) == 1):
                fg_transects_todo[j] = array([fg_transects_todo[j]])
            if (set(list(fg_transects_todo[0])) != set(list(fg_transects_todo[j]))):
                print 'can not compare transects of these scenarios. '
                print 'because they do not have the same transect numbers'
                print 'fg_transects_todo[0] was: ',fg_transects_todo[0]
                print 'j and fg_transects_todo[j] were: ',j,fg_transects_todo[j]
                sys.exit(1)

        #All transects have the same numbers, so continue
        print 'the fg_transects_todo[0] was: '
        print fg_transects_todo[0]
        print ' '
        for iwhich in fg_transects_todo[0]:       #all fg_transects_todo are equal as a set
            print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            print '++++++++++++ PROCESSING Transect number (iwhich) +++++++ ',iwhich
            print ' '

            hazard_curves_file=[]; OutZetafile=[]; #transect iwhich, multiple scenarios
            fixed_grid_file=[];                    #Scenarios could be both fine and coarse
                                                   #and coarse resolution.

            for j in range(num_scenarios):
                outdir=scenariodir[j]+'/'+'transect_fg'+str(iwhich)
                hazard_curves_file.append(outdir + '/hazardcurves_probs.npy')
                OutZetafile.append(outdir + '/' + 'ZetaFloods.npy')
                print 'The hazard curves .npy file for jth scenario when j= ',j,' was: '
                print hazard_curves_file[j]
                print ' '
                print 'The Flood zetapbar .npy file for jth scenario when j= ',j,' was: '
                print OutZetafile[j]
                print ' '
                type_name=scenarioresolution[j]
                fg_file = runsdir+ type_name +\
                                       '_xyB0_fg'+str(iwhich) +'.npy'
                fixed_grid_file.append(fg_file)
                print 'The fixed_grid_file .npy file for the jth scenario when j= ',j,' was: '
                print fixed_grid_file[j]
                print ' '

            #The plotter will automatically choose evenly distributed points for
            #the hazard curves.  If these are good enough, choose whichrows=None;
            #otherwise, make whichrows a list (say of gauge indices) and the plotter
            #will append the evenly distributed points to the whichrows list.
            #BUT NOTE THAT the gauge indices can be transect dependent, and if this
            #is a priority, we would read them in from right scenario directory,
            #and make_scenario.py would have to be modified to put them there.

            whichrows=None
            depth=True; speed=False; mflux=False;

            #CHANGE zeta and whichzetas if desired
            zeta = list(linspace(0,12,1201)) + list(linspace(12.2,24,60))
            zeta = array(zeta)
            ### corresponding to zeta=0.0,1.0,2.0,3.0,4.0,6.0,8.0,10.0,12.0
            whichzetas=[0,100,200,300,400,600,800,1000,1200]

            whichzetas=array(whichzetas)
            pbarlist=[.1,.2,.3,.4,.5,.6,.7,.8,.9,.95]   #10 zeta-contours or 10 zeta-transects

            comparisondir=projectdir + '/comparisons_'+old_new+'/comparison_'+\
                          str(comparison_number) + '/transect_fg' + str(iwhich) 
            print 'The comparison outputs for this transect will be written at: '
            print comparisondir
            print ' '
            print 'CALLING hazard_curves_transect_compare for fg transect: ',iwhich
            print ' '

            if 1:
                hazard_curves_transect_compare(hazard_curves_file,fixed_grid_file,\
                      OutZetafile,\
                      zeta,whichzetas,comparisondir,pbarlist,\
                      CCtopo.X,CCtopo.Y,CCtopo.Z,scenarionames,scenarioresolution,\
                      whichrows=whichrows,speed=speed,depth=depth,mflux=mflux)

    if (itype == 2):              #Do all the 2Dgrid type of fixed grids
        found=1
        fg_2Dgrids_todo=[]        #A list of arrays
        for j in range(num_scenarios):
            fg_2Dgrids_todo.append(loadtxt(fg_2Dgrids_todo_file[j],dtype=integer))
            if (size(fg_2Dgrids_todo[j]) == 1):
                fg_2Dgrids_todo[j] = array([fg_2Dgrids_todo[j]])
            if (set(list(fg_2Dgrids_todo[0])) != set(list(fg_2Dgrids_todo[j]))):
                print 'can not compare 2Dgrids of these scenarios. '
                print 'because they do not have the same 2Dgrid numbers'
                print 'fg_2Dgrids_todo[0] was: ',fg_2Dgrids_todo[0]
                print 'j and fg_2Dgrids_todo[j] were: ',j,fg_2Dgrids_todo[j]
                sys.exit(1)

        #All 2Dgrids have the same numbers, so continue
        print 'the fg_2Dgrids_todo[0] was: '
        print fg_2Dgrids_todo[0]
        print ' '

        for iwhich in fg_2Dgrids_todo[0]:       #all fg_2Dgrids_todo are equal as a set
            print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            print '++++++++++++ PROCESSING 2Dgrid number (iwhich) +++++++ ',iwhich
            print ' '

            hazard_curves_file=[]; OutZetafile=[]; #2Dgrid iwhich, multiple scenarios
            fixed_grid_file=[]
            for j in range(num_scenarios):
                outdir=scenariodir[j]+'/'+'2Dgrid_fg'+str(iwhich)
                hazard_curves_file.append(outdir + '/hazardcurves_probs.npy')
                OutZetafile.append(outdir + '/' + 'ZetaFloods.npy')
                print 'The hazard curves .npy file for jth scenario when j= ',j,' was: '
                print hazard_curves_file[j]
                print ' '
                print 'The Flood zetapbar .npy file for jth scenario when j= ',j,' was: '
                print OutZetafile[j]
                print ' '
                type_name=scenarioresolution[j]
                fg_file = runsdir+ type_name +\
                                       '_xyB0_fg'+str(iwhich) +'.npy'
                fixed_grid_file.append(fg_file)
                print 'The fixed_grid_file .npy file for the jth scenario when j= ',j,' was: '
                print fixed_grid_file[j]
                print ' '

            ###Note, need to read Npts1 and Npts2 from some file region_dimensions.txt
            ###      in outdir that make_scenerio.py put there.
            ###      The file has length two with the values  Npts1 and  Npts2
            ###
            ###All the scenarios will have the same Npts1 and Npts2, so just retrieve
            ###this information from the first one.

            outdir=scenariodir[0]+'/'+'2Dgrid_fg'+str(iwhich)
            dimension_file=outdir + '/' + 'region_dimensions.txt'
            dimensions=loadtxt(dimension_file,dtype=integer)
            Npts1=dimensions[0]; Npts2=dimensions[1];
            print 'The 2Dgrid was Npts1 x Npts2 points in longxlat dimensions'
            print Npts1,Npts2
            print ' '

            #The plotter will automatically choose evenly distributed points for
            #the hazard curves.  If these are good enough, choose whichrows=None;
            #otherwise, make whichrows a list (say of gauge indices) and the plotter
            #will append the evenly distributed points to the whichrows list.
            #BUT NOTE THAT the gauge indices can be 2Dgrid dependent, and if this
            #is a priority, we would read them in from right scenario directory,
            #and make_scenario.py would have to be modified to put them there.

            whichrows=None
            depth=True; speed=False; mflux=False;

            #CHANGE zeta and whichzetas if desired
            zeta = list(linspace(0,12,1201)) + list(linspace(12.2,24,60))
            zeta = array(zeta)

            ## corresponding to zeta=0.0,1.0,2.0,3,0,4.0,6.0,8.0,10.0,12.0
            whichzetas=[0,100,200,300,400,600,800,1000,1200] 

            whichzetas=array(whichzetas)
            pbarlist=[.1,.2,.3,.4,.5,.6,.7,.8,.9,.95]   #10 zeta-contours 

            comparisondir=projectdir + '/comparisons_'+old_new+'/comparison_'+\
                          str(comparison_number) + '/2Dgrid_fg' + str(iwhich) 
            print 'The comparison outputs for this 2Dgrid will be written at: '
            print comparisondir
            print ' '
            print 'CALLING hazard_curves_2Dgrid_compare for fg 2Dgrid: ',iwhich
            print ' '

            if 1:
                hazard_curves_2Dgrid_compare(hazard_curves_file,fixed_grid_file,\
                      OutZetafile,\
                      zeta,whichzetas,comparisondir,pbarlist,Npts1,Npts2,\
                      CCtopo.X,CCtopo.Y,CCtopo.Z,scenarionames,scenarioresolution,\
                      whichrows=whichrows,speed=speed,depth=depth,mflux=mflux)

if (transect_found == 0):
    print 'No transects were found to compare'

if (found == 0):
    print 'No 2Dgrids were found to compare'

