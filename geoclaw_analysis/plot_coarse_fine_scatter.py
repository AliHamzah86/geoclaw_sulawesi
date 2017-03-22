"""
Plot scatter plots comparing various things.

Need to first run compute_etamax.py 
to create files with etamin, etamean, etamax
"""

from pylab import *
import os

save_png = True
if save_png:
    png_files = os.path.abspath('../png_files')
    os.system('mkdir -p %s' % png_files)

all_runs_dir = os.path.abspath('../geoclaw_output/all_runs_npy_files')

rdict = {}

Mwlist = [8.6,8.8,9.0,9.2]
#Mwlist = [8.6,8.8,9.0]
Mw_colors = {8.6:'m', 8.8:'g', 9.0:'b', 9.2:'r'}

for Mw in Mwlist:
    for resolution in ['fine', 'coarse']:
        for fgno in [1,2,3]:
            fname = '%s_%s_etamax_fg%s.txt' % (resolution, Mw, fgno)
            d = loadtxt(os.path.join(all_runs_dir,fname), skiprows=1)
            rdict[(Mw, resolution, fgno, 'etamin')] = d[:,1]
            rdict[(Mw, resolution, fgno, 'etamean')] = d[:,2]
            rdict[(Mw, resolution, fgno, 'etamax')] = d[:,3]
            rdict[(Mw, resolution, fgno, 'hmax_sum')] = d[:,4]
            

for fgno in [1,2,3]:

    figure(500+fgno)
    clf()
    qoi = 'etamean'
    #qoi = 'hmax_sum'
    qoi_min = inf
    qoi_max = -inf

    for Mw in Mwlist:
        qoi_c = rdict[(Mw,'coarse',fgno,qoi)]
        qoi_f = rdict[(Mw,'fine',fgno,qoi)]
        #qoi_c = where(qoi_c>1., qoi_c, 1.)
        #qoi_f = where(qoi_f>1., qoi_f, 1.)
        #loglog(qoi_c,qoi_f,'o',color=Mw_colors[Mw], label='Mw = %s' % Mw)
        plot(qoi_c,qoi_f,'o',color=Mw_colors[Mw], label='Mw = %s' % Mw)
        qoi_max = max(qoi_max, nanmax(qoi_c), nanmax(qoi_f))
        qoi_min = min(qoi_min, nanmin(qoi_c), nanmin(qoi_f))
        #print "+++ fgno = %s, Mw = %s, qoi = %s" % (fgno,Mw,qoi)
        #print "+++ qoi_min = %g" % qoi_min
        #print "+++ qoi_max = %g" % qoi_max

    axis('scaled')
    #if qoi_min < 0.1: qoi_min = -0.5
    qoi_min = -1
    qoilim = [0.9*qoi_min,1.1*qoi_max]
    #loglog(qoilim, qoilim, 'k')
    plot(qoilim, qoilim, 'k')
    xlim(qoilim)
    ylim(qoilim)

    legend(loc='lower right')
    #plot([4,14],[4,14],'k')
    xlabel('computed from coarse run')
    ylabel('computed from fine run')
    title('%s over onshore points from fgmax grid %s' % (qoi,fgno))


    if save_png:
        fname = os.path.join(png_files, 'scatter_%s_fg%s.png' % (qoi,fgno))
        savefig(fname)
        print "Created ",fname

# scatter plots of etamax vs. etamean for all runs:

figureno = 600
for resolution in ['coarse','fine']:
    for fgno in [1,2,3]:
        figureno += 1
        figure(figureno)
        clf()
        qoi_min = inf
        qoi_max = -inf
        for Mw in Mwlist:
            etamean = rdict[(Mw, resolution, fgno, 'etamean')]
            etamax = rdict[(Mw, resolution, fgno, 'etamax')]
            plot(etamean,etamax,'o',color=Mw_colors[Mw], label='Mw = %s' % Mw)
            qoi_max = max(qoi_max, nanmax(etamean), nanmax(etamax))
            qoi_min = min(qoi_min, nanmin(etamean), nanmin(etamax))
        axis('scaled')
        qoi_min = -1
        qoilim = [0.9*qoi_min,1.1*qoi_max]
        plot(qoilim, qoilim, 'k')
        xlim(qoilim)
        ylim(qoilim)

        legend(loc='lower right')
        plot([0,qoi_max],[0,qoi_max],'k')
        xlabel('etamean on %s grid' % resolution)
        ylabel('etamax on %s grid' % resolution)
        title('Max vs mean eta (%s) from fgmax grid %s' % (resolution,fgno))

