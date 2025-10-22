"""

N O R T H   S U L A W E S I   T S U N A M I  
configuration : fgmax & fgout
===========================================

Plot fgmax output from GeoClaw run.

"""

import matplotlib.pyplot as plt
import os
from numpy import where, nan, linspace, nanmin, nanmax, floor
from clawpack.geoclaw import fgmax_tools
from clawpack.visclaw import geoplot

def plot_fgmax_grid(colorbar=None, zeta_lim=None, arrival_time=True, verbose=False):

    outdir = '_output'

    fg = fgmax_tools.FGmaxGrid()
    fg.outdir = outdir
    data_file = os.path.join(outdir, 'fgmax_grids.data')
    fg.read_fgmax_grids_data(fgno=1, data_file=data_file)
    fg.read_output()

    zeta = where(fg.B>0, fg.h, fg.h+fg.B)   # surface elevation in ocean
    zeta = where(fg.h>0, zeta, nan)  # mask land
    zeta_max = nanmax(zeta)
    zeta_min = nanmin(zeta)

    zeta_lim = floor(zeta_max) if zeta_lim is None else zeta_lim
    zeta_lim = max(zeta_lim, 2.0)  # ensure a minimum limit
    clines_zeta = [0.0, 0.01] + list(linspace(0.05,0.3,6)) \
                + [0.5,1.0, zeta_lim]

    plt.figure(1)
    plt.clf()

    # plot coastline:
    plt.contour(fg.X,fg.Y,fg.B,[0.],colors='k')  # coastline
    
    # plot surface elevation:
    colors = geoplot.discrete_cmap_1(clines_zeta)
    contourf = plt.contourf(fg.X,fg.Y,zeta,
                            levels=clines_zeta,
                            colors=colors,
                            extend='max', )
    if colorbar:
        cbar = plt.colorbar(contourf, label='m/s')
        tick_values = clines_zeta
        cbar.set_ticks(tick_values)
        cbar.set_ticklabels([f"{v:.2f}" for v in tick_values])
    else:      
        # Set aspect ratio to square
        plt.gca().set_aspect('equal', adjustable='box')  
    
    # plot arrival time contours and label:
    if arrival_time:
        arrival_t = fg.arrival_time/3600.  # arrival time in hours
        clines_t = linspace(0,8,17)  # hours
        clines_t_label = clines_t[2::2]  # which ones to label 
        clines_t_colors = ([.5,.5,.5],)
        con_t = plt.contour(fg.X,fg.Y,arrival_t, 
                            clines_t,
                            colors=clines_t_colors, 
                            linewidths=2.,) 
        plt.clabel(con_t, clines_t_label, colors=clines_t_colors)

    # fix axes:
    plt.ticklabel_format(style='plain',useOffset=False)
    xticks = linspace(fg.X.min(), fg.X.max(), 5)
    yticks = linspace(fg.Y.min(), fg.Y.max(), 5)
    plt.xticks(xticks, labels=[str(int(round(x))) for x in xticks], rotation=90)
    plt.yticks(yticks, labels=[str(int(round(y))) for y in yticks])
    plt.tick_params(axis='both', direction='in', top=True, right=True)

    fname = "fgmax_time_travel.png"
    plt.savefig(fname, dpi=300)
    print("Created ",fname)

    if verbose:
        print(f"maximum zeta = {zeta_max:.3f} m")
        print(f"minimum zeta = {zeta_min:.3f} m")
        print(f"last contour level:{clines_zeta[-1]:.3f} m")
        print(f"last arrival time contours: {clines_t[-1]:.1f} hours")

if __name__=="__main__":

    zeta_lim = 2.0
    plot_fgmax_grid(colorbar=True, 
                    zeta_lim=zeta_lim, 
                    arrival_time=True, 
                    verbose=True)
