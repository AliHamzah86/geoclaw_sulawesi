"""
Compute the dtopo proxies based on the dtopo files for each run.

Uses the GeoClaw output for each run.
"""

from pylab import *
import os,sys

from clawpack.geoclaw import fgmax_tools, topotools

resolution = 'coarse'  # should not matter, we only use the dtopo files!

# Which cases to process:
Mw_list = ['8.6', '8.8','9.0','9.2']
run_list = range(100)

if 1:
    # for testing:
    Mw_list = ['9.0']
    run_list = range(2)
    
# Crescent City location:
xcc = -124.1838
ycc = 41.7456


geoclaw_output = os.path.abspath('../geoclaw_output')
os.chdir(geoclaw_output)
all_runs_dir = os.path.abspath('all_runs_npy_files')


for Mw in Mw_list:
    Mw_dir = os.path.abspath('%s_%s' % (resolution,Mw))

    print "Using dtopo files from from ",Mw_dir

    topo = topotools.Topography()
    topo.read('../DataFiles/etopo1_-126_-123_40_45_1min.asc',3)

    def load_dtopo(rundir):
        from clawpack.geoclaw.dtopotools import DTopography
        dtopo = DTopography()
        dtopo.read(os.path.join(Mw_dir,rundir,'dtopo.tt3'), dtopo_type=3)
        return dtopo

    dtopo = load_dtopo('run_0')
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

    def PotentialEnergy(dZr):
        dy = 1./60. * 111.e3  # m
        dx = dy * cos(topo.Y * pi/180.)  # m
        grav = 9.81  # m/s^2
        rho_water = 1000  # kg/m^3
        eta = ma.masked_where(topo.Z>0, dZr)
        Energy = sum(eta**2 * dx * dy) * grav * rho_water * 1e-15  # PetaJoules
        return Energy

    def DisplacedVolume(dZr):
        dy = 1./60. * 111.e3  # m
        dx = dy * cos(topo.Y * pi/180.)  # m
        eta = ma.masked_where(topo.Z>0, dZr)
        volume = sum(eta * dx * dy) * 1e-9  # km^3
        return volume

    def dtopo_max_offshore(dZr):
        eta = ma.masked_where(topo.Z>0, dZr)
        return eta.max()

    fname = os.path.join(all_runs_dir, 'dtopo_proxies_%s.txt' % Mw)
    proxyfile = open(fname,'w')
    proxyfile.write('  run       dzCC (m)      eta_max (m)      PE (pJ)      dVolume\n')


    for runno in run_list:

        rundir = 'run_%s' % str(runno)

        dtopo = load_dtopo(rundir)
        dzCC = dZ_CrescentCity(dtopo)
        dZr = dtopo.dZ[0,:,:]
        PE = PotentialEnergy(dZr)
        volume = DisplacedVolume(dZr)
        eta_max_offshore = ma.masked_where(topo.Z>0, dZr).max()

        proxyfile.write('%4i   %13.8f  %13.8f  %13.8f  %13.8f\n' \
              % (runno, dzCC, eta_max_offshore, PE, volume))


    proxyfile.close()
    print "Created ",os.path.abspath(fname)
