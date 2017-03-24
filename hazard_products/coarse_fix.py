#### This program creates the coarse_mod .npy files from the coarse .npy files
#### and the coarse and fine bathymetry files.  An output text file called
#### output.coarse_fix_all has comparison information.
####
#### The user chooses below which coarse .npy files are being used as input,
#### and the directory where the output coarse_mod files are to be stored (directories
#### for all the runs have been created first with program make_coarse_mod_dirs.py).
#### If coarse_mod files are to be made for all Mw=8.6, 8.8, 9.0, and 9.2, the program
#### will need to be run 4 times, choosing the proper files below.
####
#### The associated fine .npy files are used for comparison purposes only if desired.
####
#### To run:  (Example to make the coarse_mod_9.2_runs)
####   In directory $FEMA/redclaw_1106/all_runs_npy_files/coarse_mod_9.2_runs:
####
####              python $FEMA/ptha_tools/coarse_fix.py > output.coarse_fix_all 
####
####                     (coarse_fix.py stored in $FEMA/ptha_tools)
from numpy import *
import os,sys
projectdir = os.environ['FEMA']

################   USER INPUT ###################################
#runsdir = 'all_runs_npy_files_new'
runsdir = 'redclaw_1106/all_runs_npy_files'

typedir = 'coarse_9.2_runs'
finetypedir = 'fine_9.2_runs'
out_typedir = 'coarse_mod_9.2_runs'

if 1:  #The run numbers for which coarse_mod files are to be made
    num_runs = 100
    run_nums=[]
    for j in range(num_runs):
        run_nums.append(j)

if 0:
    num_runs=1
    run_nums=[87]

if 1: #If comparison to fine is to be done
    testing=1

if 0: #No comparison to fine is to be done
    testing=0
################## END OF USER INPUT ###############################

def fixer(hmax,BathyC,BathyF):
    ####
    #### Note: If BathyF=BathyC, hmax_mod remains as hmax from below
    ####
    ####  First find etaThresh.  This will be the largest eta value where hmax and BathyC >0
    ####  (hmax + BathyC) 

    etaThresh=amax(where(logical_and(hmax > 0.0,BathyC>0.0),hmax+BathyC,0.0))
    DB=BathyF-BathyC
    hmax_mod=where(logical_and(BathyF>BathyC,hmax>DB),hmax-DB,hmax)
    hmax_mod=where(logical_and(BathyF>BathyC,hmax<=DB),0.0,hmax_mod)
    hmax_mod=where(logical_and(BathyF<BathyC,hmax>0.0),hmax-DB,hmax_mod)

    #### Prepare the fix to use for the BathyF<BathyC and hmax==0.0 case
    #### There is hmax=0 at BathyC, but a bathymetry below BathyC might also
    #### give hmax=0.0.  We are looking for hmax_mod + BathyF = etaThresh.
    #### If BathyF >= etaThresh set hmax_mod==0 (leaves hmax==0 alone).
    #### If BathyF < etaThresh set hmax_mod=etaThresh-BathyF  (when etaThresh <=BathyC)
    ####                           hmax_mod=BathyC-BathyF     (when etaThresh > BathyC)
    #### Note: when etaThresh > BathyC, can't add beyond BathyC, last where ensures this

    threshcondition = logical_and(BathyF<etaThresh,etaThresh<=BathyC)
    thresh2 =         logical_and(BathyF<BathyC,etaThresh>BathyC)
    hmax_mod=where(logical_and(threshcondition,hmax==0.0),etaThresh-BathyF,hmax_mod)
    hmax_mod=where(logical_and(thresh2,hmax==0.0),BathyC-BathyF,hmax_mod)
    return hmax_mod,etaThresh


for fg in [1,2,3]:
    dc=load(projectdir + '/' + runsdir + '/coarse_xyB0_fg' + str(fg) +'.npy')
    Bc[fg]=dc[:,2]
    df=load(projectdir + '/' + runsdir + '/fine_xyB0_fg' + str(fg) +'.npy')
    Bf[fg]=df[:,2]
#####


print 'Ready to start the modification of the coarse runs'
print ' ' 

etaThresh={}
etaThresh[1]=[]; etaThresh[2]=[]; etaThresh[3]=[];
for jj in range(num_runs):
    print '#######'
    j=run_nums[jj]
    print 'FOR run number: ',j
    filedir = projectdir + '/' + runsdir + '/' + typedir + '/run_' + str(j)

    if (testing == 1):
        filefinedir = projectdir + '/' + runsdir + '/' + finetypedir + '/run_' + str(j)

    ### Now march over the 3 fixed grids in this run and create the modified
    ### hmax file and save in dictionary etaThresh, each entry a separate list.
    for fg in [1,2,3]:
        hmax=load(filedir + '/' + 'hmax_fg' + str(fg) + '.npy')    #old coarse file
        BathyC = Bc[fg]; BathyF = Bf[fg];
        hmax_mod,etaT = fixer(hmax,BathyC,BathyF)         #hmax_mod is new coarse file
        etaThresh[fg].append(etaT)
        outdir = projectdir + '/' + runsdir + '/' + out_typedir +\
                 '/run_' + str(j) + '/' + 'hmax_fg' + str(fg) + '.npy'
        save(outdir,hmax_mod)
        print ' '
        print 'FOR fg= ',fg
        print 'The fg= ',fg,' coarse hmax file for run_ ',j,' has been modified and stored'
        print 'in: ',outdir
        print ' '

        if (testing == 1):
            hfine=load(filefinedir + '/' + 'hmax_fg' + str(fg) + '.npy')  #fine file
            diff_fine_minus_coarse = hfine - hmax
            diff_fine_minus_coarse_mod = hfine - hmax_mod
            print 'The max(abs(hfine-hmax)) -- fine to coarse: ',max(abs(hfine-hmax))
            print 'The mean(abs(hfine-hmax)) -- fine to coarse: ',mean(abs(hfine-hmax))
            print 'The max(abs(hfine-hmax_mod)) -- fine to modified coarse: ',\
                   max(abs(hfine-hmax_mod))
            print 'The mean(abs(hfine-hmax_mod)) -- fine to modified coarse: ',\
                   mean(abs(hfine-hmax_mod))
            print ' '
            diff_tocoarse_land=max(where(BathyF >= 0.0,abs(diff_fine_minus_coarse),0.0))
            diff_tocoarsemod_land=max(where(BathyF >= 0.0,abs(diff_fine_minus_coarse_mod),0.0))
            print 'The max fine to coarse on land: ',diff_tocoarse_land
            print 'The max fine to coarsemod on land: ',diff_tocoarsemod_land
            diff_tocoarse_land=mean(where(BathyF >= 0.0,abs(diff_fine_minus_coarse),0.0))
            diff_tocoarsemod_land=mean(where(BathyF >= 0.0,abs(diff_fine_minus_coarse_mod),0.0))
            print 'The mean fine to coarse on land: ',diff_tocoarse_land
            print 'The mean fine to coarsemod on land: ',diff_tocoarsemod_land
            print ' '

#### Print some things: 
print 'etaThresh is the largest eta (hmax+BathyC) where BathyC and  hmax(coarse) were > 0.0'
print ' '
print 'run_num etaThresh(fg=1) etaThresh(fg=2) etaThresh(fg=3)'
print ' '
for jj in range(num_runs):
    j=run_nums[jj]
    print j,etaThresh[1][jj],etaThresh[2][jj],etaThresh[3][jj]
    print ' '
print ' ' 
