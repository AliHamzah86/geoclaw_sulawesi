"""
Compute etamax, etamean etc.  needed by plot_comparisons.py
"""

from pylab import *
import os,sys

from clawpack.geoclaw import fgmax_tools

# Which cases to process:
resolution_list = ['coarse','fine']
Mw_list = ['8.6', '8.8','9.0','9.2']
run_list = range(100)

if 1:
    # for testing:
    resolution_list = ['coarse']
    Mw_list = ['9.0']
    run_list = range(2)

all_runs_dir = os.path.abspath('../geoclaw_output/all_runs_npy_files')
B0dir = all_runs_dir

for Mw in Mw_list:
    for resolution in resolution_list:
        topdir = '%s/%s_%s_runs' % (all_runs_dir,resolution,Mw)
        print "Processing data from %s" % topdir

        # Crescent City location:
        xcc = -124.1838
        ycc = 41.7456


        for fgno in [1,2,3]:
            etamax_fname = '%s/%s_%s_etamax_fg%s.txt' \
                        % (all_runs_dir,resolution,Mw,fgno)
            etamax_file = open(etamax_fname,'w')
            etamax_file.write('  run     min eta    mean eta    max eta (on shore)  hmax sum\n')

            fname = os.path.join(B0dir, '%s_xyB0_fg%s.npy' % (resolution,fgno))
            xyB0 = load(fname)
            X = xyB0[:,0]
            Y = xyB0[:,1]
            B0 = xyB0[:,2]

            for runno in run_list:
                fname = os.path.join(topdir, 'run_%s/hmax_fg%s.npy' \
                        % (runno,fgno))
                try:
                    hmax = load(fname)
                except:
                    print "*** Missing: %s %s run %s" % (Mw,resolution,runno)
                    hmax = emin = emean = emax = hmax_sum = nan

                if hmax is not nan:
                    # define eta as surface where wet, nan elsewhere:
                    eta = where(hmax > 0, B0+hmax, nan)
        
                    # define eta_shore as eta where wet and onshore, nan elsewhere
                    eta_shore = where(logical_and(hmax>0, B0>0), eta, nan)
        
                    # define hmax_sum as sum of hmax values on shore only
                    # (note, not multiplied by area or length)
                    hmax_sum = sum(where(B0>0, hmax, 0.))
        
                    #print "loaded %s" % fname
        
                    #print "Run %s: max eta on shore is %7.3f" \
                    #   % (runno, nanmax(eta_shore))
                    
                    emin = nanmin(eta_shore)
                    emax = nanmax(eta_shore)
                    emean = nanmean(eta_shore)
                    if isnan(emin): emin = 0.  
                    if isnan(emax): emax = 0.
                    if isnan(emean): emean = 0.

                    #print "*** Found: %s %s run %s" % (Mw,resolution,runno)
                    #import pdb; pdb.set_trace()
                    

                etamax_file.write('%4i    %8.3f    %8.3f    %8.3f    %17.3f\n' \
                        % (runno, emin, emean, emax, hmax_sum))

            etamax_file.close()
            print "Created ",etamax_fname
