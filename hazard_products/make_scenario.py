#####  To begin, a project directory (called projectdir) already exists with a subdirectory  ####
#####  called scenarios_1106 where this program can add one or more scenario directories.    ####
#####  The master weights for the KL probabilities have been copied from a previous step to  ####
#####  the file called projectdir + '/scenario_prb_wgts.txt'.  If master uniform probs are   ####
#####  to be used instead of KL probabilities, they must have been stored already in file    ####
#####  projectdir + '/scenario_prb_wgts_uniform.txt'.
#####                                                                                        ####
#####  This pgm creates a directory specified by scenarioname if it doesn't already exist    ####
#####  within the directory projectdir + '/scenarios_1106'.  Here, the projectdir is $FEMA   ####
#####  which has been set by a local environment variable using the os command.              ####
#####  Within directory scenarioname, up to 3 subdirecories are created.  Here they are      ####
#####  called 2Dgrid_fg3, transect_fg1, and transect_fg2 and the results of the main         ####
#####  probability programs are to be stored in each of these directories.                   ####
#####  The scenario probability weights are stored in the scenarioname                       ####
#####  directory after being calculated by this program, or set manually into this           ####
#####  program (events_probs) under USER INPUT for some of the centroid strategies. They are ####
#####  computed if given=False using master weights above, or set manually if given=True.    ####
#####  If given=True, the master weights file is not read in.                                ####
#####   
#####  When given=True, the lists master_run_numbers, master_magnitudes, and events_probs    ####
#####  are set manually. To get these values, the program cluster_event_probs_make.py must   ####
#####  have been already run (which uses the clusters.txt file previously calculated) to     ####
#####  find the values in these lists.  Then these lists are copied below manually, taking   ####
#####  care that the magnitude*10 is copied as just magnitude.                               #### 
#####
#####  The USER has to always specify the run numbers and their magnitudes below to make the ####
#####  scenario.  Also the event_probs are either set or computed for each scenario. As an   ####
#####  example, the program below created the directory scenario_all_SVD_40, with only 1     #### 
#####  subdirectory for the 2Dgrid_fg3 results. The master_weights_file used the KL          ####
#####  probabilities which are in file (projectdir + '/scenario_prb_wgts.txt'), and was for  ####
#####  all 400 realizations, 100 of each of the magnitudes 8.6, 8.8, 9.0, and 9.2.           ####
#####
#####  Once, the scenarioname directory is created by this program, the next steps           ####
#####  in the process will be to execute program hc_make_haz_run.py to create the            ####
#####  hazard curves .npy file for each of the subdirectories specified and then             ####
#####  execute program hc_plot_run.py to create the plots for this scenario and to           ####
#####  create the ZetaFloods.npy file in the subdirectories that is used by the              ####
#####  program hc_compare_run.py that compares various scenarios.                            ####
#####                                                                                        ####
#####  To execute:  python make_scenario.py                                                  ####
#####                                                                                        ####

import os,sys
from numpy import array, sum, zeros, linspace, loadtxt, savetxt
from numpy import sort, argsort, shape
from numpy import logical_and
projectdir = os.environ['FEMA']
commentsfile = projectdir + '/' + 'output.make_scenario_temp'
outfile = open(commentsfile,'w')

#############  USER INPUT: SET THE NEW SCENARIO PARAMETERS HERE   ###############
given=False          #Choose false if the event_probs are computed here and not manually set here.

loaded_master=0      #Have not yet loaded the master probability file. Will do so
                     #later in the program if given remains False.

### Set whether old or new run and the mag_factors used for this scenario
if 0:  #This was used when we only had 100 not 400 realizations.
    old_new='old'
    mag_factors={}
    mag_factors[9.0]=1.0;

if 1: #Choose this flag for the results for the FEMA 1106 report based on the redclaw_1106 data
    old_new='new'
    mag_factors={}
    mag_factors[9.2]=.1;  mag_factors[9.0]=.3;
    mag_factors[8.8]=.3; mag_factors[8.6]=.3;

################################   NEW 400 fine and 400 coarse data set runs ######
###
### Choose the scenarioname for this job run.  Only one such if statement should say if 1:
### All the scenarioname choices used during this project are left below as a history of
### how the scenarios were set up.
###
if 0:
    scenarioname = 'scenario_all_fine'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 1:
    #scenarioname = 'scenario_all_pseudo'
    #scenarioname = 'scenario_all_pseudo_40'
    #scenarioname = 'scenario_all_pseudoCM'
    #scenarioname = 'scenario_all_pseudoCM_40'
    #scenarioname = 'scenario_all_pseudoCM_dzCC_PE'
    #scenarioname = 'scenario_all_pseudoCM_dzCC_PE_40'
    scenarioname  = 'scenario_all_SVD_40'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0:
    scenarioname = 'scenario_all_pseudo_dtopo_addto'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0: 
    #scenarioname = 'scenario_all_pseudo_dtopo'
    scenarioname = 'scenario_all_pseudo_dtopo_40'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0: 
    scenarioname = 'scenario_all_pseudo_dtopo_uniform'
    master_weights_file = projectdir + '/scenario_prb_wgts_uniform.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0: 
    scenarioname = 'scenario_all_pseudo_dtopo_uniform_addto'
    master_weights_file = projectdir + '/scenario_prb_wgts_uniform.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0:
    scenarioname = 'scenario_all_fine_uniform'
    master_weights_file = projectdir + '/scenario_prb_wgts_uniform.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0:
    scenarioname = 'scenario_all_pseudo_uniform'
    master_weights_file = projectdir + '/scenario_prb_wgts_uniform.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0:
    scenarioname = 'scenario_all_coarse'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0:
    scenarioname = 'scenario_all_coarse_mod'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in range(100):
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0:
    scenarioname = 'scenario_highest10_each'
    highest_numbers=[43,19,25,26,58,83,13,20,59,99]
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]; master_magnitudes=[];
    for mag in [8.6,8.8,9.0,9.2]:
        for irun in highest_numbers:
            master_run_numbers.append(irun)
            master_magnitudes.append(mag)
if 0: 
    if 1:
        scenarioname = 'scenario_highest20_uniform'
        master_weights_file = projectdir + '/scenario_prb_wgts_uniform.txt'
        iend=-21
    if 0:
        scenarioname = 'scenario_highest40_overall'
        master_weights_file = projectdir + '/scenario_prb_wgts.txt'
        iend=-41
    if 0:
        scenarioname = 'scenario_highest20_overall'
        master_weights_file = projectdir + '/scenario_prb_wgts.txt'
        iend=-21
    if 0:
        scenarioname = 'scenario_highest40_uniform'
        master_weights_file = projectdir + '/scenario_prb_wgts_uniform.txt'
        iend=-41
    #
    #  Determine the highest 20 or highest 40
    master_probs = loadtxt(master_weights_file)
    loaded_master=1
    all_weights=zeros(400); all_indices=range(400);
    all_weights[0:100]=master_probs*mag_factors[8.6]
    all_indices[0:100]=range(100);                      #run numbers
    all_weights[100:200]=master_probs*mag_factors[8.8]
    all_indices[100:200]=range(100);                    #run numbers
    all_weights[200:300]=master_probs*mag_factors[9.0]
    all_indices[200:300]=range(100);                    #run numbers
    all_weights[300:400]=master_probs*mag_factors[9.2]
    all_indices[300:400]=range(100);                    #run numbers
    all_indices=array(all_indices)
    all_weights_sorted=sort(all_weights)
    indices_sorted=argsort(all_weights)

    ### Now the 20 biggest or 40 biggest are on the bottom with their indices
    ### iend=-21 if the 20 biggest, iend=-41 if the 40 biggest are desired
    indiceskeep=indices_sorted[-1:iend:-1]                     
    master_run_numbers=all_indices[indiceskeep]      
    master_run_numbers=list(master_run_numbers)
    master_magnitudes=[]
    for k in indiceskeep:
        if (0 <= k <=99):
            master_magnitudes.append(8.6)
        elif (100 <= k <=199):
            master_magnitudes.append(8.8)
        elif (200 <= k <=299):
            master_magnitudes.append(9.0)
        elif (300 <= k <= 399):
            master_magnitudes.append(9.2)
        else:
            print 'indiceskeep not right'
            print 'indiceskeep was: ', indiceskeep
            print 'program stopping'
            sys.exit(1)

if 0: #for 1106
    scenarioname = 'scenario_coarse_etaCM_40cls'
    given=True
    master_run_numbers=[4, 27, 69, 64, 83, 80, 78, 74, 83, 89, 80, 17, 34, 58, 90, 85, 73,\
         52, 20, 24, 43, 25, 81, 75, 23, 38, 14, 56, 9, 3, 19, 82, 20, 43, 80, 54, 85, 71, 48, 60]
    master_magnitudes=[9.0,8.6,9.0,8.6,8.8,8.6,9.2,9.2,9.0,8.8,8.8,8.8,9.0,9.0,9.0,8.8,9.2,9.0,9.0,9.2,\
    9.0,9.2,9.0,9.0,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2]
    events_probs=[0.0070955, 0.0516308, 0.0769082, 0.0849732, 0.1119294, 0.0885116, 0.0506133, 0.0058719,\
    0.1293458, 0.0347736, 0.0320252, 0.0175379, 0.0168437, 0.0211720, 0.0149277, 0.0058945,\
    0.0125007, 0.0556975, 0.0269617, 0.0139625, 0.0380737, 0.0134969, 0.0018320, 0.0081554,\
    0.0079291, 0.0023500, 0.0009169, 0.0015614, 0.0111914, 0.0069309, 0.0129286, 0.0096150,\
    0.0060613, 0.0129183, 0.0000282, 0.0016272, 0.0015384, 0.0001202, 0.0006880, 0.0028604]

if 0: #for 1106, also uses CM
    scenarioname = 'scenario_dzCC_PE_40cls'
    given=True
    master_run_numbers=[60, 11, 18, 5, 60, 17, 18, 58, 43, 82, 48, 46, 31, 48, 46, 25, 66, 43, 25, 66,\
                        34, 34, 59, 37, 84, 87, 87, 84, 27, 4, 98, 78, 50, 4, 88, 27, 4, 67, 4, 29]
    master_magnitudes=[9.2,9.2,9.2,9.2,9.0,9.2,9.0,9.2,9.2,9.0,8.8,9.0,9.2,8.6,8.8,9.0,9.2,8.6,8.8,9.0,\
                       8.6,8.8,9.2,8.6,9.0,8.8,9.0,9.2,8.8,8.6,9.2,9.0,9.2,8.8,9.2,9.2,9.0,9.2,9.2,9.2]
    events_probs=[0.0015796, 0.0000031, 0.0020670, 0.0018093, 0.0047482, 0.0056490, 0.0068782, 0.0073533,\
    0.0173837, 0.0216976, 0.0150842, 0.0742110, 0.0146432, 0.0163770, 0.0924508, 0.0582282,\
    0.0253353, 0.1350876, 0.0595284, 0.0759056, 0.1376652, 0.1093504, 0.0059979, 0.0088586,\
    0.0365374, 0.0171918, 0.0144312, 0.0103165, 0.0043728, 0.0020116, 0.0042387, 0.0041168,\
    0.0011692, 0.0032458, 0.0006408, 0.0011495, 0.0020116, 0.0000033, 0.0006700, 0.0000006]

if 0: #for 1106
    scenarioname = 'scenario_coarse_eta_40cls'
    given=True
    master_run_numbers=[4, 79, 78, 24, 25, 29, 33, 98, 66, 89, 83, 84, 17, 34,\
                        23, 58, 8, 97, 99, 20, 26, 72, 25, 81, 9, 38, 85, 32,\
                        58, 3, 82, 19, 20, 43, 80, 54, 85, 71, 48, 60]
    master_magnitudes=[8.8,8.6,8.8,8.8,8.8,9.0,8.8,9.0,8.8,8.8,9.0,9.0,8.8,9.0,9.0,9.0,\
                       8.8,9.2,9.0,9.0,9.2,9.0,9.2,9.0,9.0,9.2,9.0,\
                       9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2]
    events_probs=[0.0070955, 0.0519678, 0.0911752, 0.0723138, 0.0684502,\
                 0.1226208, 0.0251782, 0.0527334, 0.0195534, 0.0347736,\
                 0.1154940, 0.0088830, 0.0178327, 0.0163995, 0.0073638,\
                 0.0211720, 0.0162301, 0.0089175, 0.0162843, 0.0269617,\
                 0.0853942, 0.0160625, 0.0134969, 0.0018615, 0.0025582,\
                 0.0023405, 0.0047871, 0.0020020, 0.0121346, 0.0069309,\
                 0.0122606, 0.0129286, 0.0060613, 0.0129183, 0.0000282,\
                 0.0017287, 0.0014369, 0.0001202, 0.0006911, 0.0028572]

if 0: #for 1106
    scenarioname = 'scenario_dtopo_40cls'
    given=True
    master_run_numbers=[60, 77, 8, 82, 18, 43, 42, 7, 85, 28, 99, 22, 54, 43,\
                        10, 25, 25, 23, 71, 26, 19, 19, 95, 97, 76, 36,\
                        95, 84, 35, 35, 98, 32, 83, 27, 98, 78, 37, 4, 27, 4]
    master_magnitudes=[9.2,9.2,9.2,9.2,9.0,9.2,9.2,9.2,9.2,9.2,9.0,9.2,9.0,9.0,9.0,\
                       9.2,9.0,8.8,8.6,9.0,8.8,9.2,9.0,8.6,9.2,8.6,\
                       9.2,8.8,8.8,9.2,9.0,9.2,9.2,8.8,9.2,9.0,9.2,9.0,9.2,9.2]
    events_probs=[0.0015827, 0.0009291, 0.0029471, 0.0025770, 0.0162115, 0.0188491,\
                 0.0000055, 0.0066853, 0.0014901, 0.0009973, 0.0214552, 0.0033491,\
                 0.0092026, 0.0710284, 0.0075327, 0.0082373, 0.0262047, 0.0833968,\
                 0.1814100, 0.0793379, 0.0706602, 0.0174293, 0.0024768, 0.1421500,\
                 0.0081175, 0.0535304, 0.0004406, 0.0855832, 0.0371289, 0.0000128,\
                 0.0094455, 0.0035364, 0.0065947, 0.0076320, 0.0020155, 0.0053510,\
                 0.0006308, 0.0020116, 0.0011529, 0.0006705]
if 0: #for 1106
    scenarioname = 'scenario_c_eta_20cls_uniform'
    given=True
    master_run_numbers=[4,69,12,0,45,80,49,10,90,26,24,20,6,32,9,33,20,43,77,18] 
    master_magnitudes=[8.8,8.8,8.6,8.8,8.8,8.8,8.8,9.0,9.0,9.2,9.2,9.0,9.0,9.2,9.2,9.2,\
                       9.2,9.2,9.2,9.2]
    events_probs=[0.013,0.123,0.145,0.151,0.146,0.048,0.061,0.046,0.021,0.058,0.044,0.035,\
                  0.030,0.005,0.024,0.025,0.007,0.007,0.004,0.007]

if 0: #for 1106
    scenarioname = 'scenario_coarse_eta_20cls'
    given=True
    master_run_numbers=[4,69,12,0,45,80,49,10,90,26,24,20,6,32,9,33,20,43,77,18] 
    master_magnitudes=[8.8,8.8,8.6,8.8,8.8,8.8,8.8,9.0,9.0,9.2,9.2,9.0,9.0,9.2,9.2,9.2,\
                       9.2,9.2,9.2,9.2]
    events_probs=[0.0070955,\
                  0.1132282,\
                  0.1022286,\
                  0.1984964,\
                  0.1099523,\
                  0.0667654,\
                  0.1054736,\
                  0.0259920,\
                  0.0200689,\
                  0.0956334,\
                  0.0322829,\
                  0.0367547,\
                  0.0091128,\
                  0.0020229,\
                  0.0205595,\
                  0.0284910,\
                  0.0060895,\
                  0.0129183,\
                  0.0031656,\
                  0.0036685]


if 0:
    scenarioname = 'scenario_c_eta_20cls_uniform'
    given=True
    master_run_numbers=[4,35,50,52,88,66,70,17,29,66,7,72,85,1,26,57,95,12,82,77]
    master_magnitudes=[8.6,8.6,8.6,8.6,8.8,8.8,8.8,8.8,9.0,9.0,9.0,9.0,9.0,9.2,\
                       9.2,9.2,9.2,9.2,9.2,9.2]
    events_probs=[0.018, 0.075, 0.087, 0.126, 0.114, 0.018, 0.084, 0.087, 0.087,\
                  0.018, 0.076, 0.08, 0.023, 0.018, 0.012, 0.02, 0.006, 0.024,\
                  0.015, 0.012]

if 0: #for 1106, also uses CM
    scenarioname = 'scenario_dzCC_PE_20cls'
    given=True
    master_run_numbers=[60, 8, 11, 53, 49, 46, 25, 43, 46, 26, 34, 13, 34, 84, 69, 20, 69, 2, 4, 4]
    master_magnitudes=[9.2,9.2,9.0,9.2,9.0,9.0,9.2,8.8,8.6,9.0,8.6,9.2,8.8,9.0,8.8,9.2,9.0,9.2,9.0,9.2]
    events_probs=[0.0019795, 0.0030486, 0.0158511, 0.0300500, 0.0298662,\
                  0.0742110, 0.0194094, 0.1398382, 0.1278694, 0.1199356,\
                  0.1701090, 0.0265670, 0.1407766, 0.0646864, 0.0073626,\
                  0.0155639, 0.0073626, 0.0016779, 0.0031645, 0.0006705]
if 0:
    scenarioname = 'scenario_coarse_etaCM_20cls'
    given=True
    master_run_numbers=[29,50,49,44,83,86,80,10,52,82,20,49,77,32,17,33,20,43,77,18]
    master_magnitudes=[8.6,8.6,8.6,8.8,9.0,9.2,8.8,9.0,9.0,9.0,9.0,9.0,9.0,9.2,9.2,9.2,9.2,9.2,9.2,9.2]
    events_probs=[0.0598751, 0.127646, 0.2022356, 0.0780631, 0.1397848, 0.0638356, 0.034471,\
                  0.0257666, 0.0569326, 0.0389451, 0.0388129, 0.0518749, 0.0084597, 0.0020229,\
                  0.02168, 0.0237522, 0.0060895, 0.0129183, 0.0034568, 0.0033773]

if 0:
    scenarioname = 'scenario_coarse_eta_20cls'
    given=True
    master_run_numbers=[4,35,50,52,17,66,70,88,7,29,66,72,85,1,\
                        12,26,57,77,82,95]
    master_magnitudes=[8.6,8.6,8.6,8.6,8.8,8.8,8.8,8.8,9.0,9.0,\
                       9.0,9.0,9.0,9.2,9.2,9.2,9.2,9.2,9.2,9.2]
    events_probs=[0.00384440,\
                  0.08044660,\
                  0.07047680,\
                  0.15622320,\
                  0.08200300,\
                  0.01357820,\
                  0.08083160,\
                  0.11494920,\
                  0.05158593,\
                  0.09831640,\
                  0.02261080,\
                  0.11078627,\
                  0.00672627,\
                  0.01627293,\
                  0.03050840,\
                  0.00780007,\
                  0.02293327,\
                  0.00630787,\
                  0.01905393,\
                  0.00474487]

if 0:
    scenarioname = 'scenario_dzCC_dtopo_20cls_uniform'
    given=True
    master_run_numbers=[60,43,48,62,43,52,76,34,71,3,43,83,3,66,57,63,86,2,27,4]
    master_magnitudes=[9.2,9.2,9.2,9.2,9.0,9.0,9.2,9.2,8.8,\
                       9.0,8.8,9.2,8.6,9.0,9.2,9.2,8.8,8.6,9.0,9.2]
    events_probs=[0.006,0.008, 0.011, 0.016, 0.039, 0.075, 0.01, 0.019, 0.105, 0.081,\
                 0.093, 0.008, 0.189, 0.054, 0.009, 0.01, 0.105, 0.102, 0.057, 0.003]

if 0: #for 1106, uniform
    scenarioname = 'scenario_dzCC_dtopo_20cls_uniform'
    given=True
    master_run_numbers=[18,77,18,43,92,62,43,26,39,7,86,61,80,80,59,98,37,27,27,4]
    master_magnitudes=[9.2,9.2,9.0,9.2,9.2,9.2,9.0,9.2,9.0,8.8,9.2,8.6,8.8,9.0,9.2,9.2,9.2,\
                       9.0,9.2,9.2]
    events_probs=[0.006,0.003,0.02,0.021,0.002,0.017,0.066,0.01,0.07,0.147,0.026,0.318,0.15,0.089,\
                 0.019,0.003,0.004,0.024,0.003,0.002]

if 0: #for 1106
    scenarioname = 'scenario_dzCC_dtopo_etamax_20cls'
    given=True
    master_run_numbers=[18,77,18,43,92,62,43,26,39,7,86,61,80,80,59,98,37,27,27,4]
    master_magnitudes=[9.2,9.2,9.0,9.2,9.2,9.2,9.0,9.2,9.0,8.8,9.2,8.6,8.8,9.0,9.2,9.2,9.2,\
                       9.0,9.2,9.2]
    events_probs=[0.0045299,\
                  0.0009291,\
                  0.0162115,\
                  0.0255344,\
                  0.0000557,\
                  0.0060915,\
                  0.0848064,\
                  0.0204734,\
                  0.0498478,\
                  0.1538788,\
                  0.0199175,\
                  0.3217412,\
                  0.1530632,\
                  0.1160525,\
                  0.0146697,\
                  0.0020283,\
                  0.0007137,\
                  0.0076320,\
                  0.0011529,\
                  0.0006705]

if 0:
    scenarioname = 'scenario_dzCC_dtopo_etamax_20cls'
    given=True
    master_run_numbers=[2,3,43,71,86,3,27,43,52,66,4,34,43,48,\
                        57,60,62,63,76,83]
    master_magnitudes=[8.6,8.6,8.8,8.8,8.8,9.0,9.0,9.0,9.0,9.0,\
                       9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2,9.2]
                        
    events_probs=[0.10924140,\
                  0.19264180,\
                  0.11442360,\
                  0.06739200,\
                  0.11579700,\
                  0.09594780,\
                  0.02583740,\
                  0.07513720,\
                  0.02555300,\
                  0.07802880,\
                  0.00108193,\
                  0.01893007,\
                  0.01856027,\
                  0.00241740,\
                  0.00446540,\
                  0.00325620,\
                  0.01408420,\
                  0.00700753,\
                  0.02551360,\
                  0.00468340]

####################################  OLD BELOW HERE #############################
### All scenarios below this line were used in testing or with the old 100 data set
### Set scenarioname, the master indices of the run numbers used for this scenario
#scenarios_old/scenario_highest10
if 0:
    scenarioname = 'scenario_highest10'
    master_run_numbers=[43,19,25,26,58,83,13,20,59,99]
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
if 0:
    scenarioname = 'scenario_highest10_newzetas'
    master_run_numbers=[43,19,25,26,58,83,13,20,59,99]
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'

#scenarios_old/scenario_uniform10
if 0:
    scenarioname = 'scenario_uniform10'
    master_weights_file = projectdir + '/scenario_prb_wgts_uniform.txt'
    master_run_numbers=[43,19,25,26,58,83,13,20,59,99]

#scenarios_old/scenario_all_fine_uniform
if 0:
    scenarioname = 'scenario_all_fine_uniform'
    master_weights_file = projectdir + '/scenario_prb_wgts_uniform.txt'
    master_run_numbers=[]
    for irun in range(100):
        master_run_numbers.append(irun)

#scenarios_old/scenario_all_fine
if 0:
    scenarioname = 'scenario_all_fine'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]
    for irun in range(100):
        master_run_numbers.append(irun)

#scenarios_old/scenario_all_fine_newzetas
if 0:
    scenarioname = 'scenario_all_fine_newzetas'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]
    for irun in range(100):
        master_run_numbers.append(irun)

#scenarios_old/scenario_all_coarse
if 0:
    scenarioname = 'scenario_all_coarse'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]
    for irun in range(100):
        master_run_numbers.append(irun)

#scenarios_old/scenario_all_coarse_newzetas
if 0:
    scenarioname = 'scenario_all_coarse_newzetas'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]
    for irun in range(100):
        master_run_numbers.append(irun)

#scenarios_old/scenario_all_coarse_mod
if 0:
    scenarioname = 'scenario_all_coarse_mod'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]
    for irun in range(100):
        master_run_numbers.append(irun)

#scenarios_old/scenario_h10_coarse_mod
if 0:
    scenarioname = 'scenario_h10_coarse_mod'
    master_weights_file = projectdir + '/scenario_prb_wgts.txt'
    master_run_numbers=[]
    for irun in range(100):
        master_run_numbers.append(irun)

### Set the fixed grid transects and fixed grid  2Dgrid numbers to do for this scenario
if 0:                                            #Doing both transects and 2D grids
    fg_types_todo=[1,2]                          #1 means 1D transects, 2 means 2D grids
    fg_transects_todo=[1,2]                      #master fgnos
    fg_2Dgrids_todo=[3]                          #master fgnos
    list1=[88,78]
    fg_2Ddimensions={}
    fg_2Ddimensions[3]=list1                     #Npts1,Npts2  (Npts1 in long, Npts2 in lat)
                                                 #A dictionary of lists. Dictionary is keyed
                                                 #by the fixed grid number. Will just have
                                                 #keys corresponding to the 2D fixed grids.
   
if 1:                                            #Doing just 2D grid, for DonSub's SVD
    fg_types_todo=[2]                            #1 means transects, 2 means 2D grid
    fg_2Dgrids_todo=[3]                          #master fgnos: 1 and 2 were transects, 3 was 2D grid
    list1=[88,78]
    fg_2Ddimensions={}
    fg_2Ddimensions[3]=list1

#####
if 0:
    scenariodir = projectdir + '/scenarios_' + old_new + '/' + scenarioname

if 1: #Always choose for the FEMA 1106 results
    scenariodir = projectdir + '/scenarios_1106' + '/' + scenarioname
####################################################################
########################  END OF USER INPUT ########################

### FOR the old_new == 'old' runs, the runs_todo_file only has one column (the run numbers)
### FOR the old_new == 'new' runs, the runs_todo_file has two columns.
### The  1st colm is run number, 2nd is magnitude directory                          
### For example, run_10 with magnitude 8.8 would have  10  88  in the file           
### The 88 would locate the right directory to pull run 10 from.                  

runs_todo_file = scenariodir + '/runs_todo.txt'

scenario_weights_file=scenariodir + '/scenario_prb_wgts.txt'
fg_types_todo_file = scenariodir + '/fg_types_todo.txt'
fg_transects_todo_file = scenariodir + '/fg_transects_todo.txt'
fg_2Dgrids_todo_file = scenariodir + '/fg_2Dgrids_todo.txt'
#
commentsfile2 = scenariodir + '/' + 'output.make_scenario'

outfile.write(' The project directory is: %s ' %projectdir)
outfile.write(' \n')
outfile.write(' The scenario name is: %s ' %scenarioname)
outfile.write(' \n')
outfile.write(' \n')
outfile.write(' The scenario directory desired is: %s ' %scenariodir)
outfile.write(' \n')
outfile.write(' \n')
outfile.write(' The conditional weights for this scenario are in file: \n')
outfile.write(' %s ' %scenario_weights_file)
outfile.write(' \n')
outfile.write(' \n')

### Create the scenario directory if it doesn't exist
os.system('mkdir -p %s' % scenariodir)
outfile.write(' The desired scenario directory was created called: \n')
outfile.write(' %s ' %scenariodir)
outfile.write(' \n')
outfile.write(' \n')

### Write the master indices of the run numbers used for this scenario
### in the scenariodir on the runs_todo.txt file
outfile.write(' The run numbers in master_probs to use in this scenario were: \n')
outfile.write(' %s ' % master_run_numbers)
outfile.write(' \n')
outfile.write(' \n')

if (old_new == 'old'):
    run_mag=zeros((len(master_run_numbers),1))
    run_mag[:,0]=master_run_numbers
    outfile.write(' The run numbers were stored on file: \n')
elif (old_new == 'new'):
    run_mag=zeros((len(master_run_numbers),2))
    run_mag[:,0]=master_run_numbers
    run_mag[:,1]=array(master_magnitudes)*10
    outfile.write(' The run numbers (and magnitudes) were stored on file: \n')
else:
    print 'old_new was not set properly, program exiting'
    sys.exit(1)

savetxt(runs_todo_file,run_mag,fmt='%4i')
#savetxt(runs_todo_file,run_mag)
outfile.write(' %s ' %runs_todo_file)
outfile.write(' \n')
outfile.write(' \n')

if (old_new == 'old'):
    outfile.write(' The magnitudes of the master_run numbers were all 9.0: \n')
    outfile.write(' \n')
else:
    outfile.write(' The magnitudes of the master_run numbers were: \n')
    outfile.write(' %s ' % master_magnitudes)
    outfile.write(' \n')
    outfile.write(' \n')

outfile.write(' The mag_factors dictionary was: \n')
outfile.write(' %s ' % mag_factors)
outfile.write(' This dictionary not used when given==True, given was: %s ' %given)
outfile.write(' \n')
outfile.write(' \n')

### Compute the probabilities to use in this scenerio and store them
### in the scenariodir directory in the scenario_weights_file 
###
if (logical_and(loaded_master != 1,given==False)):
    master_probs = loadtxt(master_weights_file)
##


#### New code, Uses different weights when have different magnitudes        #####
if (logical_and(old_new == 'new',given==False)):  #calculate events_probs if not given
    sum=0.0; events_probs=[];
    for k in range(len(master_run_numbers)):
        index=master_run_numbers[k]
        sum += master_probs[index]*mag_factors[master_magnitudes[k]]
    for k in range(len(master_run_numbers)):
        index=master_run_numbers[k]
        events_probs.append(master_probs[index]*mag_factors[master_magnitudes[k]]/sum)

#### Old code, only had one magnitude which was 9.0                         #####
if (old_new == 'old'):
    sum=0.0; events_probs=[];
    for index in master_run_numbers:  
        sum += master_probs[index]
    for index in master_run_numbers:
        events_probs.append(master_probs[index]/sum)

events_probs=array(events_probs)
allprobs=events_probs.sum()
outfile.write(' The sum of all the events_probs was: %s ' %allprobs)
outfile.write(' \n') 
outfile.write(' \n') 

if (old_new == 'old'):
    outfile.write(' The run number and conditional probs for this scenerio are: \n')
    j=0
    for index in master_run_numbers:
        outfile.write('%s  %s ' %(index, events_probs[j]))
        outfile.write(' \n')
        j += 1
if (old_new == 'new'):
    outfile.write(' The run number, magnitude, and cond probs for this scenerio are: \n')
    j=0
    for index in master_run_numbers:
        outfile.write('%s  %s %s' %(index, master_magnitudes[j], events_probs[j]))
        outfile.write(' \n')
        j += 1

outfile.write(' \n')
outfile.write(' \n')
savetxt(scenario_weights_file,events_probs)
outfile.write(' These probabilities were stored in file: \n')
outfile.write(' %s ' %scenario_weights_file)
outfile.write(' \n')
outfile.write(' \n')

### Loop over the fg_types_todo and create an appropriate output directory
### in the scenario directory for each fg_transect and each fg_2Dgrid. Also
### save the fg_transects_todo_file, the fg_2Dgrids_todo_file, and the
### fg_types_todo file.
###
for itype in fg_types_todo:
    if (itype == 1):               #set up directories for all the fg_transects
        for iwhich in fg_transects_todo:
            tdir = scenariodir + '/' + 'transect_fg' + str(iwhich)  
            os.system('mkdir -p %s' % tdir)
            outfile.write(' made directory: %s ' %tdir)
            outfile.write(' \n')
        savetxt(fg_transects_todo_file,fg_transects_todo,fmt='%2i')
    if (itype == 2):               #set up directories for all the fg_2Dgrids
        for iwhich in fg_2Dgrids_todo:
            tdir = scenariodir + '/' + '2Dgrid_fg' + str(iwhich)  
            os.system('mkdir -p %s' % tdir)
            outfile.write(' made directory: %s ' %tdir)
            outfile.write(' \n')
            dimension_file=tdir + '/' + 'region_dimensions.txt'
            dimensions = fg_2Ddimensions[iwhich] 
            outfile.write(' The 2Dgrid was Npts1 x Npts2 -long,lat- \n')
            outfile.write(' %s %s ' %(dimensions[0],dimensions[1]))
            outfile.write(' \n')
            outfile.write(' \n')
            savetxt(dimension_file,dimensions,fmt='%5i')
        savetxt(fg_2Dgrids_todo_file,fg_2Dgrids_todo,fmt='%2i')
    outfile.write(' \n')
    outfile.write(' \n')
savetxt(fg_types_todo_file,fg_types_todo,fmt='%1i')
outfile.close()

### Move the commentsfile to the scenario directory
os.system('mv  %s %s' %(commentsfile,commentsfile2)) 
print 'Job run output stored at: '
print commentsfile2 
print ' '
