"""
Plot the fgmax grid results

"""

import matplotlib
matplotlib.use('Agg')  # for non-interactive

from pylab import *
import os,sys

from clawpack.geoclaw import fgmax_tools
from clawpack.geoclaw.topotools import Topography
from clawpack.geoclaw.geoplot import discrete_cmap_1

save_png = True
if save_png:
    png_files = os.path.abspath('../png_files')
    os.system('mkdir -p %s' % png_files)


all_runs_dir = os.path.abspath('../geoclaw_output/all_runs_npy_files')


# Topography for contour lines:
CCtopo = Topography()
CCtopo.read('../DataFiles/CCtopo_1_3sec.tt3',3)

# Crescent City location:
xcc = -124.1838
ycc = 41.7456

fgno = 3

cres = {'fine': 'r', 'coarse': 'b'}  # color to use for each resolution


# Which cases to process:
resolution_list = ['coarse','fine']
Mw_list = ['8.6', '8.8','9.0','9.2']
run_list = range(100)

if 0:
    # for testing:
    resolution_list = ['coarse']
    Mw_list = ['9.0']
    run_list = range(2)

figure(1000,figsize=(10,8))
clf()

for Mw in Mw_list:
    for runno in run_list:

        #figure(1000+runno,figsize=(10,8)) # to open all figs
        clf()

        gray = '#888888'
        contour(CCtopo.X, CCtopo.Y, CCtopo.Z, arange(0,21,2), colors=gray)

        #for resolution in ['fine','coarse']:
        for resolution in resolution_list:
            topdir = '%s/%s_%s_runs' % (all_runs_dir,resolution,Mw)
            #fname = os.path.join(topdir, '%s_xyB0_fg%s.npy' % (resolution,fgno))
            fname = '%s/%s_xyB0_fg3.npy' % (all_runs_dir,resolution)
            xyB0 = load(fname)
            X = xyB0[:,0]
            Y = xyB0[:,1]
            B0 = xyB0[:,2]
            X = reshape(X, (88,78), order='F')
            Y = reshape(Y, (88,78), order='F')
            B0 = reshape(B0, (88,78), order='F')

            if False and resolution == 'fine':
                contour(X, Y, B0, arange(0,12.5,.5), colors='g')  # to compare

            fname = os.path.join(topdir, 'run_%s/hmax_fg%s.npy' \
                    % (runno,fgno))
            hmax = load(fname)
            Hmax = reshape(hmax, (88,78), order='F')
            Hmax = where(B0>0, Hmax, nan)
            #contourf(X, Y, Hmax, [0,5,20], colors=['b','r'], alpha=0.5)
            #contourf(X, Y, Hmax, [0.1,2,4,100], colors=['b','#ffaaaa','r'], alpha=0.8)

            clines_h = [.01] + list(arange(0.5,10.5,.5))
            colors_h = discrete_cmap_1(clines_h)
            contourf(X, Y, Hmax, clines_h, colors=colors_h, alpha=0.95, extend='max')
            colorbar()

            if 0:
                contour(X, Y, Hmax, [0.01], colors=cres[resolution], \
                        linewidths=3, label=resolution)
                contour(X, Y, Hmax, [5.], colors=cres[resolution], \
                        linewidths=1, label=resolution)

            # define eta as surface where wet, nan elsewhere:
            eta = where(Hmax > 0, B0+Hmax, nan)

            eta_shore = where(logical_and(Hmax>0, B0>0), eta, nan)
            # define eta_shore as eta where wet and onshore, nan elsewhere

            #print "loaded %s" % fname

            print "Run %s: max eta on shore is %7.3f" \
               % (runno, nanmax(eta_shore))
            print "     max h on shore is %7.3f" % nanmax(Hmax)
            
        legend(loc='lower right',fontsize=10)
        gca().set_aspect(1./cos(41.7*pi/180.))
        ticklabel_format(format='plain',useOffset=False)
        xticks(rotation=20)
        x1 = X.min(); x2 = X.max(); y1 = Y.min(); y2 = Y.max()
        plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'g-')
        #axis([x1,x2,y1,y2])  # to crop to fgmax region alone

        title('Run %s Mw = %s' % (runno,Mw))

        if save_png:
            png_dir = '%s/%s_%s/run_%s/' \
                    % (png_files,resolution, Mw, runno)
            os.system('mkdir -p %s' % png_dir)
            fname = 'run%s_Mw%s_fg%s.png' % (runno, Mw, fgno)
            fname = os.path.join(png_dir, fname)
            savefig(fname, bbox_inches='tight')
            print "Created ",fname
