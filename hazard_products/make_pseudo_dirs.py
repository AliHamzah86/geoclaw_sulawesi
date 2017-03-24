##### Makes all the pseudo run directories for a particular pseudo strategy
##### User needs to choose below which pseudoname is being used for this run.
#####
##### To run:   python make_pseudo_dirs.py  (run from $FEMA/ptha_tools directory)
#####
import os,sys
projectdir = os.environ['FEMA']

#####   USER INPUT, CHOOSE pseudoname  #######################
if 0:
    pseudoname = 'pseudo_40_'
if 0:
    pseudoname = 'pseudo_dtopo_40_'

########  Below are ones using coarse mod (CM) to make the pseudo
########  All these below are for the KL probabilities. Those without
########  the 40 extension have 20 clusters.
if 0:
    pseudoname = 'pseudoCM_'
if 0:
    pseudoname = 'pseudoCM_40_'
if 0:
    pseudoname = 'pseudoCM_dzCC_PE_'
if 1:
    pseudoname = 'pseudoCM_dzCC_PE_40_'
#####
#####   END OF USER INPUT #############

#dir86 = projectdir + '/all_runs_npy_files_new/' + pseudoname + '8.6_runs'
#dir88 = projectdir + '/all_runs_npy_files_new/' + pseudoname + '8.8_runs'
#dir90 = projectdir + '/all_runs_npy_files_new/' + pseudoname + '9.0_runs'
#dir92 = projectdir + '/all_runs_npy_files_new/' + pseudoname + '9.2_runs'

dir86 = projectdir + '/redclaw_1106/all_runs_npy_files/' + pseudoname + '8.6_runs'
dir88 = projectdir + '/redclaw_1106/all_runs_npy_files/' + pseudoname + '8.8_runs'
dir90 = projectdir + '/redclaw_1106/all_runs_npy_files/' + pseudoname + '9.0_runs'
dir92 = projectdir + '/redclaw_1106/all_runs_npy_files/' + pseudoname + '9.2_runs'

for irun in range(100):
    tdir =  dir86 + '/run_' +str(irun) 
    os.system('mkdir -p %s' % tdir)
    print ' made directory: %s ' %tdir
    tdir =  dir88 + '/run_' +str(irun) 
    os.system('mkdir -p %s' % tdir)
    print ' made directory: %s ' %tdir
    tdir =  dir90 + '/run_' +str(irun) 
    os.system('mkdir -p %s' % tdir)
    print ' made directory: %s ' %tdir
    tdir =  dir92 + '/run_' +str(irun) 
    os.system('mkdir -p %s' % tdir)
    print ' made directory: %s ' %tdir
    print ' '
