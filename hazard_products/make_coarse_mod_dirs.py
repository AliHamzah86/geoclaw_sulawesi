##### Makes all the coarse_mod directories
#####
##### To run:    python make_coarse_mod_dirs.py  (run from $FEMA/ptha_tools)
#####
import os,sys
projectdir = os.environ['FEMA']

#dir86 = projectdir + '/all_runs_npy_files_new/coarse_mod_8.6_runs'
#dir88 = projectdir + '/all_runs_npy_files_new/coarse_mod_8.8_runs'
#dir90 = projectdir + '/all_runs_npy_files_new/coarse_mod_9.0_runs'
#dir92 = projectdir + '/all_runs_npy_files_new/coarse_mod_9.2_runs'

dir86 = projectdir + '/redclaw_1106/all_runs_npy_files/coarse_mod_8.6_runs'
dir88 = projectdir + '/redclaw_1106/all_runs_npy_files/coarse_mod_8.8_runs'
dir90 = projectdir + '/redclaw_1106/all_runs_npy_files/coarse_mod_9.0_runs'
dir92 = projectdir + '/redclaw_1106/all_runs_npy_files/coarse_mod_9.2_runs'

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
