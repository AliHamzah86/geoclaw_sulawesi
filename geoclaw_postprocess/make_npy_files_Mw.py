
"""
Run this code to make fine_runs or coarse_runs npy files
from the fgmax and gauge files stored from each run.  

"""

from pylab import *
import os,sys,glob

from clawpack.geoclaw import fgmax_tools
import clawpack.pyclaw.gauges as gauges
from scipy.interpolate import interp1d

# Which cases to process:
resolution_list = ['coarse','fine']
Mw_list = ['8.6', '8.8','9.0','9.2']
run_list = range(100)

if 1:
    # for testing:
    resolution_list = ['coarse']
    Mw_list = ['9.0']
    run_list = range(2)

make_gauge_npy_files = False

geoclaw_output = os.path.abspath('../geoclaw_output')

os.chdir(geoclaw_output)
print "Working in directory: ", geoclaw_output

os.system('mkdir -p all_runs_dir')

all_runs_dir = os.path.abspath('all_runs_npy_files')

# Crescent City location:
xcc = -124.1838
ycc = 41.7456


for resolution in resolution_list:
    for Mw in Mw_list:

        Mw_dir = '%s_%s' % (resolution,Mw)

        def load_dtopo(rundir):
            from clawpack.geoclaw.dtopotools import DTopography
            dtopo = DTopography()
            dtopo.read(os.path.join(rundir,'dtopo.tt3'), dtopo_type=3)
            return dtopo

        dtopo = load_dtopo('%s/run_0' % Mw_dir)
        i1cc = find(dtopo.x<xcc).max()
        j1cc = find(dtopo.y<ycc).max()
        a1cc = (xcc-dtopo.x[i1cc])/(dtopo.x[i1cc+1]-dtopo.x[i1cc])
        a2cc = (ycc-dtopo.y[j1cc])/(dtopo.y[j1cc+1]-dtopo.y[j1cc])
        if (a1cc<0.) or (a1cc>1.) or (a2cc<0.) or (a2cc>1.):
            print '*** Interpolation to CC not correct!'

        def dZ_CrescentCity(dtopo):
            dZr = dtopo.dZ[0,:,:]
            dzy1 = (1.-a1cc)*dZr[j1cc,i1cc] + a1cc*dZr[j1cc,i1cc+1]
            dzy2 = (1.-a1cc)*dZr[j1cc+1,i1cc] + a1cc*dZr[j1cc+1,i1cc+1]
            dzcc = (1.-a2cc)*dzy2 + a2cc*dzy1
            return dzcc



        new_runs_dir = os.path.join(all_runs_dir, '%s_%s_runs'  % (resolution,Mw))
        os.system('mkdir -p %s' % new_runs_dir)

        dzfile = open(os.path.join(new_runs_dir, 'dzCC.txt'), 'w')


        for runno in run_list:

            rundir = 'run_%s' % str(runno)
            fullrundir = os.path.join(Mw_dir, rundir)
            outdir = os.path.join(fullrundir,'_output')
            if not os.path.isfile(os.path.join(fullrundir,'_output/fort.FG1.aux1')):
                print "Missing files in %s" % rundir
                continue

            newrundir = '%s/%s' % (new_runs_dir, rundir)
            os.system('mkdir -p %s' % newrundir)
            print "Created ",newrundir

            dtopo = load_dtopo(fullrundir)
            dzCC = dZ_CrescentCity(dtopo)
            dzfile.write('%4i   %13.8f\n' % (runno, dzCC))
            fname = os.path.join(newrundir, 'dzCC.txt')
            savetxt(fname, array([dzCC]))
            print "Created ",fname

            for fgno in [1,2,3]:
                fg = fgmax_tools.FGmaxGrid()
                fg.read_input_data('fgmax%s_%s.txt' % (fgno, resolution))
                fg.read_output(fgno,os.path.join(fullrundir,'_output'))
                hmax = reshape(fg.h, -1, order='F')  # as 1d array
                fname = os.path.join(newrundir, 'hmax_fg%s.npy' % fgno)
                save(fname, array(hmax))
                print "Created ",fname

            if make_gauge_npy_files:
                gauge_files = glob.glob(fullrundir + '/_output/gauge0*txt')
                gaugenos = []
                for f in gauge_files:
                    gaugenos.append(int(f[-9:-4]))
                if runno==0: print "found gauges ",gaugenos

                tt = linspace(0,8900,891)  # equal times for sampling gauges
                for gaugeno in gaugenos:
                    g = gauges.GaugeSolution(gaugeno, outdir)
                    eta_func = interp1d(g.t,g.q[3,:],bounds_error=False)
                    eta_t = eta_func(tt)
                    h_func = interp1d(g.t,g.q[0,:],bounds_error=False)
                    h_t = h_func(tt)
                    h_eta = vstack((h_t,eta_t)).T
                    fname = os.path.join(newrundir, 'gauge%s.npy' % str(gaugeno).zfill(5))
                    save(fname,h_eta)
                    #print "Created ",fname
                print "Created gauge npy files"

                    
        dzfile.close()
