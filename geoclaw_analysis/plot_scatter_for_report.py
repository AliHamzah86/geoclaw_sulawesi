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

cell_area = 30.*cos(41.7*pi/180.)  # area of fgmax grid cell in m^2

for Mw in Mwlist:
    for resolution in ['fine', 'coarse']:
        for fgno in [1,2,3]:
            fname = '%s_%s_etamax_fg%s.txt' % (resolution, Mw, fgno)
            d = loadtxt(os.path.join(all_runs_dir,fname), skiprows=1)
            rdict[(Mw, resolution, fgno, 'etamin')] = d[:,1]
            rdict[(Mw, resolution, fgno, 'etamean')] = d[:,2]
            rdict[(Mw, resolution, fgno, 'etamax')] = d[:,3]
            rdict[(Mw, resolution, fgno, 'hmax_sum')] = d[:,4]
            rdict[(Mw, resolution, fgno, 'volume')] = d[:,4] * cell_area
            

fgno = 3

#qoi_list = [('etamean',[-0.5,10,-0.5,10])]
#qoi_list = [('etamax',[-0.5,10,-0.5,10])]
#qoi_list = [('volume',[-10e3,300e3,-10e3,300e3])]

qoi_list = [('etamean',[-0.5,10,-0.5,10]), ('etamax',[-0.5,10,-0.5,10]), 
            ('volume',[-10e3,300e3,-10e3,300e3])]

for qoi_item in qoi_list:
    figure(500, figsize=(9,9))
    clf()
    qoi = qoi_item[0]
    qoi_min = inf
    qoi_max = -inf

    for Mw in Mwlist:
        qoi_c = rdict[(Mw,'coarse',fgno,qoi)]
        qoi_f = rdict[(Mw,'fine',fgno,qoi)]
        plot(qoi_c,qoi_f,'o',color=Mw_colors[Mw], label='Mw = %s' % Mw)
        qoi_max = max(qoi_max, nanmax(qoi_c), nanmax(qoi_f))
        qoi_min = min(qoi_min, nanmin(qoi_c), nanmin(qoi_f))

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
    xlabel('computed from coarse run',fontsize=20)
    ylabel('computed from fine run',fontsize=20)
    title('%s from fgmax grid %s' % (qoi,fgno), fontsize=20)

    #ticklabel_format(style='sci',scilimits=(-1,10))
    if qoi == 'volume':
        xticks(rotation=20)

    if save_png:
        fname = 'scatter_%s_fg%s.png' % (qoi,fgno)
        fname = os.path.join(png_files, fname)
        savefig(fname)
        print "Created ",fname

    for limits in qoi_item[1:]:
        axis(limits)

    #ticklabel_format(style='sci',scilimits=(-1,10))

    if save_png:
        fname = 'scatter_%s_fg%s_zoom.png' % (qoi,fgno)
        fname = os.path.join(png_files, fname)
        savefig(fname)
        print "Created ",fname
        


if 1:
    # scatter plots of etamax vs. etamean for fine grid runs:

    figureno = 600
    for resolution in ['coarse','fine']:
        for fgno in [3]:
            figureno += 1
            figure(figureno, figsize=(9,9))
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
            #xlim(qoilim)
            #ylim(qoilim)
            axis([-1,40,-1,40])

            legend(loc='lower right')
            plot([0,qoi_max],[0,qoi_max],'k')
            xlabel('etamean on %s grid' % resolution, fontsize=20)
            ylabel('etamax on %s grid' % resolution, fontsize=20)
            title('Max vs mean eta (%s) from fgmax grid %s' \
                    % (resolution,fgno),fontsize=20)

        if save_png:
            fname = 'scatter_etamean_etamax_%s_fg%s.png' \
                     % (resolution,fgno)
            fname = os.path.join(png_files, fname)
            savefig(fname)
            print "Created ",fname
