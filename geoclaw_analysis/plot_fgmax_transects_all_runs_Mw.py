"""
Plot all the transect cross sections together.

"""

from pylab import *
import os,sys

from clawpack.geoclaw import fgmax_tools

save_png = True
if save_png:
    png_files = os.path.abspath('../png_files')
    os.system('mkdir -p %s' % png_files)


all_runs_dir = os.path.abspath('../geoclaw_output/all_runs_npy_files')
B0dir = all_runs_dir


# Which cases to process:
resolution_list = ['coarse','fine']
Mw_list = [8.6,8.8,9.0,9.2]
run_list = range(50)  # since all 100 makes plots too cluttered

if 0:
    # for testing:
    resolution_list = ['coarse']
    Mw_list = [9.0]
    run_list = range(2)


eta_color = 'r'


# Crescent City location:
xcc = -124.1838
ycc = 41.7456


for Mw in Mw_list:
    for fgno in [1,2]:
        for resolution in resolution_list:
            fname = os.path.join(B0dir, '%s_xyB0_fg%s.npy' % (resolution,fgno))
            xyB0 = load(fname)
            X = xyB0[:,0]
            Y = xyB0[:,1]
            B0 = xyB0[:,2]

            figure(30+fgno,figsize=(14,5))
            clf()

            for runno in run_list:
                topdir = all_runs_dir + '/%s_%s_runs' % (resolution,Mw)
                fname = os.path.join(topdir, 'run_%s/hmax_fg%s.npy' \
                        % (runno,fgno))
                try:
                    hmax = load(fname)
                except:
                    print "*** Skipping runno %s: no hmax npy file" % runno
                    continue

                # define eta as surface where wet, nan elsewhere:
                eta = where(hmax > 0, B0+hmax, nan)

                eta_shore = where(logical_and(hmax>0, B0>0), eta, nan)
                # define eta_shore as eta where wet and onshore, nan elsewhere

                #print "loaded %s" % fname

                #print "Run %s: max eta on shore is %7.3f" \
                #   % (runno, nanmax(eta_shore))
                
                B0_neg = where(B0<0, B0, nan)
                fill_between(Y,B0_neg,0,color=[.6,.6,1])
                plot(Y,B0+hmax,eta_color)
                #plot(Y,B0+hmax,label=str(runno))
                plot(Y,B0,'g')
                ticklabel_format(format='plain',useOffset=False)

            #legend(loc='lower right',fontsize=10)
            grid(True)
            title('Mw %s events, %s eta on fgmax grid %s' % (Mw,resolution,fgno))

            if save_png:
                #png_dir = '%s/%s_%s/' \
                #        % (png_files,resolution, Mw)
                png_dir = png_files
                os.system('mkdir -p %s' % png_dir)
                fname = 'Mw%2i_fg%s_%sruns_%s.png' \
                    % ((10*Mw), fgno, len(run_list), resolution)            
                fname = os.path.join(png_dir, fname)
                savefig(fname, bbox_inches='tight')
                print "Created ",fname
