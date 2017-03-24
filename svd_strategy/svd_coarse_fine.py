from __future__ import print_function

from clawpack.geoclaw import fgmax_tools
import matplotlib as mpl
mpl.use('Agg')
label_size = 14
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 

import matplotlib.pylab as pl
import os, sys
import numpy as np

from matplotlib.ticker import FuncFormatter


r"""
we assume that GeoClaw data is contained in a subdirectory structure 

    ../geoclaw_output/
                      coarse_8.6/
                                 run_*/
                      coarse_8.8/
                                 run_*/
                      coarse_9.0/
                                 run_*/
                      coarse_9.2/
                                 run_*/

                      fine_8.6/
                               run_*/
                      fine_8.8/
                               run_*/
                      fine_9.0/
                               run_*/
                      fine_9.2/
                               run_*/

"""


def read_fgmax(j, fgno=3, mesh='fine', \
                  magnitude=9.0, \
                  geoclaw_outputdir='../geoclaw_output',\
                  driver_home='../geoclaw_driver'):
    """
    read fgmax values for a given run number j and return 
    a fgmax_tools.FGmaxGrid object.

    read_fgmax reads from the directory

        ../geoclaw_output/[grid-size]_[Mw]/run_[j]

    for example,

        ../geoclaw_output/fine_8.8/run_2

    if j == -1 it reads from run_B (empty run) and returns the bathymetry,
    for example,

        ../geoclaw_output/fine_B0
    
    """

    output_subdir = '../geoclaw_output/' + str(mesh) + '_' + str(magnitude)
    fgmax_txt = os.path.join(driver_home,'fgmax' + str(fgno) + '_fine.txt')
    
    if j == -1:
        # read in bathymetry from empty run stored in run_B
        output_subdir = '../geoclaw_output/' + str(mesh) + '_B0'
        output_dir = os.path.join(output_subdir,'_output')

    else:
        # read in bathymetry from empty run stored in run_j
        rundir = 'run_' + str(j)
        output_dir = os.path.join(output_subdir,\
                                  rundir,'_output')

    fgrun = fgmax_tools.FGmaxGrid();
    fgrun.read_input_data(fgmax_txt);

    fgrun.read_output(fgno, output_dir);
    
    return fgrun


def read_dtopos(run_list=range(100),\
                rundir='../geoclaw_output/fine_9.0'):
    r"""

    for run numbers given in kwarg run_list,
    return a matrix with dtopos flattened and stacked as columns, 

    """
    
    # get dtopo from run_0 to get dimensions
    dZ = np.loadtxt(rundir + '/run_0/dtopo.tt3', skiprows=9)
    DZ = np.empty((dZ.shape[0]*dZ.shape[1] ,100))
    
    for j in run_list:
        print('reading run: ' + str(j), end="\r")
        sys.stdout.flush()
        dZ = np.loadtxt(rundir + '/run_' + str(j) + '/dtopo.tt3', skiprows=9)
        DZ[:,j] = dZ.flatten()
    print('\n- done.')

    return DZ


def read_inundation(run_list=range(100),mesh='fine',magnitude=9.0,field='h+B'):
    r"""

    read inundation values from fgmax grids corresp to run_list 

    field = 'h+B' or 'h'

    """

    # read bathymetry from empty run
    fg_runB = read_fgmax(-1,fgno=3, mesh=mesh)
    B0 = fg_runB.B.flatten()

    dims = fg_runB.B.shape      # get dimensions
    
    N = len(run_list)
    A = np.empty((np.prod(dims),N))
    Dry_bathy = np.ones(np.prod(dims))
    
    # construct matrix of inundations and find dry portions of domains
    for j in range(N):
        print('reading: '+ str(magnitude) + ' ' + mesh + ' ' + str(j), end="\r")
        sys.stdout.flush()
        fg_runj = read_fgmax(j, fgno=3, mesh=mesh, magnitude=magnitude);
        if field == 'h+B':
            A[:,j] = (fg_runj.h + fg_runB.B).flatten()
        elif field == 'h':
            A[:,j] = (fg_runj.h).flatten()
        Dry_bathy = Dry_bathy * (fg_runj.h.flatten() < 1e-8)
    print('\n- done.')

    # get mask for dry parts of the domain
    mask = np.logical_and(B0 > 0, Dry_bathy < 1)

    return A, mask, dims


def compute_singvals():
    r"""

    compute singular values

    """

    mesh0 = 'fine'
    Mw0 = 9.0
    M = 6

    A, mask, dims = read_inundation(mesh=mesh0, magnitude=Mw0)

    # naive SVD
    A1 = A.copy()
    U1,s1,V1 = np.linalg.svd(A1)

    # normalized SVD, each column has N(0,1)
    A2 = A.copy()
    Astd = A2.std(axis=0)
    Amean = A2.mean(axis=0)
    A2 = (A2 - Amean) / A2.std(axis=0)
    
    U2,s2,V2 = np.linalg.svd(A2)

    # direct low-rank approx
    A1p = np.dot(np.dot(U1[:,:M],np.diag(s1[:M])),V1[:M,:])

    # scaled low-rank approx
    A2p = np.dot(np.dot(U2[:,:M],np.diag(s2[:M])),V2[:M,:])
    A2p = A2p*Astd + Amean

    error1 = np.linalg.norm( A.flatten() - A1p.flatten() )
    error2 = np.linalg.norm( A.flatten() - A2p.flatten() )

    print(' error for naive: ' + str(error1))
    print(' error for svd: ' + str(error2))

    return s1,s2

 
def get_sing_coords_X(Nruns = 100):
    r"""

    read in dtopo files and output its SVD

    """

    DZ = read_dtopos()
    UX,sX,VX = np.linalg.svd(DZ,full_matrices=False)
    
    return UX,sX,VX,DZ


def get_sing_coords_Y(mesh='fine',magnitude=9.0,Nruns = 100,field='h+B'):
    """
    read in fgmax data from runs and return its SVD,
    and also return the bathymetry

    """
    
    # read bathymetry
    fg_runB = read_fgmax(-1,fgno=3, mesh=mesh)
    B0 = fg_runB.B.flatten()

    dims = fg_runB.B.shape      # get dimensions
    
    A = np.empty((np.prod(dims),Nruns))
    Dry_bathy = np.ones(np.prod(dims))
    
    # construct matrix of inundations and find dry_bathy
    for j in range(Nruns):
        fg_runj = read_fgmax(j, fgno=3, mesh=mesh, magnitude=magnitude);
        if field == 'h+B':
            A[:,j] = (fg_runj.h + fg_runB.B).flatten()
        elif field == 'h':
            A[:,j] = (fg_runj.h).flatten()
        Dry_bathy = Dry_bathy * (fg_runj.h.flatten() < 1e-8)

    # zero out dry_bathy
    mask = np.logical_and(B0 > 0, Dry_bathy < 1)
     
    UY,sY,VY = np.linalg.svd(A,full_matrices=False)
    
    return UY,sY,VY,A,mask,B0


def long_lat(x, pos):
    """
    a helper function to change the tick format
    (The two args are the value and tick position)

    """
    return '%3.2f' % (x)


def read_outputs(magnitude):
    """
    read outputs from runs belonging to a certain magnitude
    """

    fg_run0 =  read_fgmax(0)

    Uc,sc,Vc,Ac,maskc,Bc = \
            get_sing_coords_Y(mesh='coarse', magnitude=magnitude);
    Uf,sf,Vf,Af,maskf,Bf = \
            get_sing_coords_Y(mesh='fine',   magnitude=magnitude);

    return fg_run0,Uc,sc,Vc,Ac,maskc,Bc,Uf,sf,Vf,Af,maskf,Bf


def compute_corr(Uf,Uc,M=7,\
        plot_fname='coarse_fine_corr.png',save_plot=False,verbose=False):
    """
    plot correlation between fine and coarse runs

    """

    cfU = np.dot(Uf.T,Uc)

    mode_error_list = []
    for k in range(M):
        mode_corr = cfU[k,k]

    if save_plot:
        pl.figure(figsize=(8,7));
        pl.pcolor(cfU[:M,:M],vmin=-1.,vmax=1.,cmap='RdBu');
        pl.title('Similarity of modes',size=16);
        pl.xlabel('coarse grid mode #', size=16);
        pl.ylabel('fine grid mode #', size=16);
        pl.colorbar();

        pl.savefig(plot_fname,dpi=100)

    return cfU


def select_fine_list(Yc,M=5,m=2):
    """
    return fine_runs with the largest *relative* components in the
    singular mode directions

    """
    fine_list = []
    Yc_s = Yc[:M,:]     # pick off M coordinates, default M = 4

    for k in range(M):
        jj = range(M)
        jj.remove(k)
        rYc_s = np.abs(Yc_s[k:(k+1),:]) / np.sum(np.abs(Yc_s[jj,:]),axis=0)
        fine_list = fine_list + \
                      np.argsort(rYc_s.flatten())[::-1][0:m].tolist()

    return fine_list


def compute_hAapprox(hUf,Yc,cfhU,hscaling,M):
    """
    compute approximation Af, hAapprox

    """

    hsgn = np.diag(np.diag(cfhU[:M,:M]))
    hsgn = hsgn*hsgn            # not using signs now - Yct handles it
    hsgn = hscaling*hsgn
    hAapprox = np.dot(hUf[:,:M],np.dot(hsgn,Yc[:M,:]))

    return hAapprox


def compute_hscaling(cfhU,Uc,hUf,Yc,Af,fine_list,maskf,lb=0.8,ub=1.2,M=4,m=60):
    """
    compute proper scaling to be used to adjust the approximate modes hUf
    do this by brute-force: for each scaling, comparing the approximation  
    with corresponding fine runs computed in fine_list

    """

    hscaling_list = np.linspace(lb,ub,m).tolist()
    error_list = []
    cfhU = compute_corr(hUf,Uc, M=M)
    print('(compute_hscaling) computing hscaling...')

    for hscaling in hscaling_list:

        hAapprox = compute_hAapprox(hUf,Yc,cfhU,hscaling,M)

        error_mean_list = []
        for k in fine_list:
            error_mean_list.append(np.mean(np.abs((hAapprox[maskf,k] - Af[maskf,k]))))

        error_list.append(np.mean(error_mean_list))
        
    ii = np.argsort(error_list)[0]
    hscaling_best = hscaling_list[ii]

    return hscaling_best


def compute_error(hAapprox,Af,Bf,maskf,fg_run0,\
                  save_plots=False, Nruns=100, S = (88,78), reltol=1e-1):
    """
    compute absolute and relative error between hAapprox

    return a list for each of the runs (default Nruns=100)
    """

    error_list = []
    error_relative_list = []

    for k in range(Nruns):

        Efk = (Af[:,k] - Bf)
        hEapprox = (hAapprox[:,k] - Bf)
    
        true_masked = np.ma.MaskedArray(\
                            Efk.reshape(S), mask=1 - maskf.reshape(S))
        approx_masked = np.ma.MaskedArray(\
                            hEapprox.reshape(S), mask=1 - maskf.reshape(S))
        error_masked = \
            np.ma.MaskedArray((Af[:,k] - hAapprox[:,k]).reshape(S),\
                                      mask=1 - maskf.reshape(S))
        error_masked_relative = np.ma.MaskedArray(\
            ((Af[:,k] - hAapprox[:,k])/Efk).reshape(S),\
                                 mask=1 - maskf.reshape(S))
        error_eta_masked = np.ma.MaskedArray(\
                Bf.reshape(S), mask=1 -  maskf.reshape(S))

        error_list.append(np.mean(np.abs(error_masked)))

        ii = (Efk > reltol) 
        error_relative_list.append(\
                np.mean(np.abs((Af[ii,k] - hAapprox[ii,k])/Efk[ii])))

    
        # plot and save errors
        if save_plots:

            ncontours = 20

            xlim0 = fg_run0.X.min() - 0.0005
            xlim1 = fg_run0.X.max() + 0.0005
            ylim0 = fg_run0.Y.min() - 0.0005
            ylim1 = fg_run0.Y.max() + 0.0005

            formatter = FuncFormatter(long_lat)
            v0 = max([Efk[maskf].max(), hEapprox[maskf].max()])
            v1 = min([Efk[maskf].min(), hEapprox[maskf].min()])

            v2 = np.ceil(10.*max([abs(v0), abs(v1)])) / 10.

            vmax = v2
            vmin = -v2

            contour_levels = np.linspace(vmin, vmax, ncontours)


            pl.figure(figsize=(14,15));

            pl.subplot(2,2,1)
            pl.contourf(fg_run0.X, fg_run0.Y,true_masked,\
                    vmin=vmin,vmax=vmax,cmap='GnBu');
            ax = pl.gca()
            ax.yaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_formatter(formatter)
            pl.xlim([xlim0, xlim1])
            pl.ylim([ylim0, ylim1])
            for label in ax.get_xticklabels():
                label.set_rotation(30)
            pl.colorbar();
            pl.title('True (fine-grid run)', size=16);

            pl.subplot(2,2,2)
            pl.contourf(fg_run0.X, fg_run0.Y, \
                        approx_masked,vmin=vmin,vmax=vmax,cmap='GnBu');
            pl.xlim([xlim0, xlim1])
            pl.ylim([ylim0, ylim1])
            ax = pl.gca()
            ax.yaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_formatter(formatter)
            pl.xlim([xlim0, xlim1])
            pl.ylim([ylim0, ylim1])
            for label in ax.get_xticklabels():
                label.set_rotation(30)
            pl.colorbar();
            pl.title('Approximation', size=16);

            pl.subplot(2,2,3)
            pl.contourf(fg_run0.X, fg_run0.Y, error_masked,\
                    contour_levels,cmap='RdBu');
            ax = pl.gca()
            ax.yaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_formatter(formatter)
            pl.xlim([xlim0, xlim1])
            pl.ylim([ylim0, ylim1])
            for label in ax.get_xticklabels():
                label.set_rotation(30)
            pl.colorbar();
            pl.title('True - Approx', size=16);

            contour_levels = np.linspace(-2.,2.,20)
            pl.subplot(2,2,4)
            pl.contourf(fg_run0.X, fg_run0.Y,error_masked_relative,\
                    contour_levels,cmap='RdBu');
            ax = pl.gca()
            ax.yaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_formatter(formatter)
            pl.xlim([xlim0, xlim1])
            pl.ylim([ylim0, ylim1])
            for label in ax.get_xticklabels():
                label.set_rotation(30)
            pl.colorbar();
            pl.title('(True - Approx) / True', size=16);
            pl.subplots_adjust(hspace=0.29,left=0.075,right=0.97,top=0.95,bottom=0.08)
            plot_fname = 'error_run_' + str(k) + '_' + str(magnitude) + '.png'
            pl.savefig(plot_fname,dpi=100)

            pl.close()

    return error_list,error_relative_list


def create_pcontour(A,pwgts,zeta0,maskf,fg_run0,\
                    S = (88,78),plot_fname='pcontour.png'):
    """
    given h + dB, and probability weights, save p-contour and values
    for given exceedance level zeta0

    """

    pcA_array = np.zeros(S)

    Nruns = A.shape[1]

    for j in range(Nruns):
        h4j = A[:,j].reshape(S)

        # add prob. if exceedance level reached.
        pcA_array += (h4j > zeta0) * pwgts[j] 

    pcA = np.ma.MaskedArray(pcA_array, mask=1 - maskf.reshape(S))

    # save plot 
    plot_title = 'p-contour, zeta = ' + str(zeta0)
    ncontours = 11
    contour_levels = np.linspace(0.,1. + 1e-4,ncontours)

    xlim0 = fg_run0.X.min() - 0.0005
    xlim1 = fg_run0.X.max() + 0.0005
    ylim0 = fg_run0.Y.min() - 0.0005
    ylim1 = fg_run0.Y.max() + 0.0005

    formatter = FuncFormatter(long_lat)

    pl.figure();
    pl.contourf(fg_run0.X, fg_run0.Y,pcA,\
                    contour_levels,cmap='GnBu');

    # format axis
    ax = pl.gca()
    ax.yaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_formatter(formatter)
    pl.xlim([xlim0, xlim1])
    pl.ylim([ylim0, ylim1])
    for label in ax.get_xticklabels():
        label.set_rotation(30)
    pl.colorbar();
    pl.title(plot_title, size=16);
    pl.savefig(plot_fname)

    pl.close()
    return pcA


def create_zetacontour(A,pwgts,p0,maskf,fg_run0,plot_fname='zetacontour.png'):
    """
    given h + dB, and probability weights, save zeta-contour and values
    for given p-level p0

    """
    S = (88,78)     # shape hard-coded in (bad)

    Npts = S[0]*S[1]
    Nruns = A.shape[1]

    zcA_array = np.zeros(Npts)
    pwgts = np.array(pwgts)     # typecast into numpy array

    for j in range(Npts):
        zeta_vals = A[j,:]
        
        zeta_sort_args = np.argsort(zeta_vals)[::-1]
        pwgts_sorted = pwgts[zeta_sort_args]
        ii = np.cumsum(pwgts_sorted) >= p0      # "left"-aligned step function
        if ii.any():
            zcA_array[j] = max([0., zeta_vals[zeta_sort_args[ii]].max()])
        else:
            # margin case, when p0 is too small, set to zero
            zcA_array[j] = 1e-12

    zcA = np.ma.MaskedArray(zcA_array.reshape(S), mask=1 - maskf.reshape(S))

    # save plot 
    plot_title = 'zeta contour, p = ' + str(p0)
    ncontours = 15

    xlim0 = fg_run0.X.min() - 0.0005
    xlim1 = fg_run0.X.max() + 0.0005
    ylim0 = fg_run0.Y.min() - 0.0005
    ylim1 = fg_run0.Y.max() + 0.0005

    formatter = FuncFormatter(long_lat)

    pl.figure();
    pl.contourf(fg_run0.X, fg_run0.Y,zcA,\
                    ncontours,cmap='GnBu');

    # format axis
    ax = pl.gca()
    ax.yaxis.set_major_formatter(formatter)
    ax.xaxis.set_major_formatter(formatter)
    pl.xlim([xlim0, xlim1])
    pl.ylim([ylim0, ylim1])
    for label in ax.get_xticklabels():
        label.set_rotation(30)
    pl.colorbar();
    pl.title(plot_title, size=16);
    pl.savefig(plot_fname)

    pl.close()
    return zcA

### compute approximate solution with SVD


## get singular mode coordinates for seafloor deformation
## UX,sX,VX,DZ = get_singular_coords_X()



def run_reconstruction(Mw,data_dict,Nruns = 100, M = 4, S = (88,78), \
                       driver_home='../geoclaw_driver',savedir='_output'):
    """
    build reconstruction for given magnitude Mw

    """

    fg_run0,Uc,sc,Vc,Ac,maskc,Bc,Uf,sf,Vf,Af,maskf,Bf = data_dict[Mw]

    prb_wgts_ffname = os.path.join(driver_home,'scenario_prb_wgts.txt')
    pwgts = np.loadtxt(prb_wgts_ffname)
    
    # compute coordinates in singular modes S*V
    Yc = np.dot(np.diag(sc),Vc)
    Yf = np.dot(np.diag(sf),Vf)

    # select fine-runs to use for approximating fine modes
    fine_list = select_fine_list(Yc,M=M)
    print(fine_list)

    data_dict['fine_list'].append(fine_list)
    
    # compute lean SVD
    hUf,hsf,hVf = np.linalg.svd(Af[:, fine_list],full_matrices=False) 
    hYf = np.dot(np.diag(hsf),hVf)
    
    # the transformation between Yc and Yf
    x = np.linalg.solve(np.dot(Yc[:M,fine_list],Yc[:M,fine_list].T), \
                        np.dot(Yc[:M,fine_list],hYf[:M,:].T))

    # plot and save the transformation matrix
    if 1:
        pl.figure()
        v0 = np.max(np.abs(x))
        pl.pcolor(x,cmap='RdBu',vmax=v0,vmin=-v0)
        pl.colorbar()

        plot_fname = os.path.join(savedir,\
                            'Y_transform_' + str(Mw) + '.png')
        pl.xticks(np.arange(M+1), [str(j) for j in range(M+1)])
        pl.yticks(np.arange(M+1), [str(j) for j in range(M+1)])
        pl.savefig(plot_fname, dpi=100)
        pl.close()

    Yct = np.dot(x.T,Yc[:M,:])

    cfhU = compute_corr(hUf,Uc,M=M)
         
    maskf = np.ones(S[0]*S[1],dtype=bool)   # temporary
    hscaling = 1.0  
    print('hscaling = ' + str(hscaling))


    hAapprox = compute_hAapprox(hUf,Yct,cfhU,hscaling,M)
    error_list, error_relative_list = compute_error(hAapprox,Af,Bf,\
                                            maskf,fg_run0,save_plots=False)
    
    print('mean error = ' + str(np.mean(error_list)))
    print('mean relative error = ' + str(np.mean(error_relative_list)))
    
    ### compute approximation by NMF (TODO)
    
    B0f = np.array(Bf.copy()) 
    B0f = B0f * (B0f >= 0.)

    hCapprox = hAapprox - (B0f).reshape((S[0]*S[1],1))
    C = Af - (B0f).reshape((S[0]*S[1],1))
    
    return hCapprox,C



def create_hazard_maps(hCapprox,C,pwgts,fg_run0,maskf,zeta0 = 1.0, p0 =0.5,\
                       map_type='zeta',savedir='_output',save_plot=True):

    if map_type=='zeta':
        print('create zeta-contours')
        # create true zeta contour
        plot_fname = os.path.join(savedir, \
                                  'zetacontour_true_' + str(p0) + '.png')
        zcA = create_zetacontour(C,pwgts,p0,maskf,fg_run0,\
                plot_fname=plot_fname)

        # create true approx zeta contour
        plot_fname = os.path.join(savedir, \
                                  'zetacontour_approx_' + str(p0) + '.png')
        zcAapprox = create_zetacontour(hCapprox,pwgts,p0,maskf,fg_run0,\
                plot_fname=plot_fname)

        # compute zeta contour error
        zcError = zcA - zcAapprox
        print('----')
        print('mean zeta-contour error (p=' + str(p0) + '): ' \
                + str(np.mean(np.abs(zcError))))
        print('max zeta-contour error (p=' + str(p0) + '): ' \
                + str(np.max(np.abs(zcError))))
        print('----')

    elif map_type=='p':
        print('create p-contours')
        plot_fname = os.path.join(savedir, \
                              'pcontour_true_' + str(zeta0) + '.png')
        pcA = create_pcontour(C,pwgts,zeta0,maskf,fg_run0,\
            plot_fname=plot_fname)
        plot_fname = os.path.join(savedir, \
                              'pcontour_approx_' + str(zeta0) + '.png')
        pcAapprox = create_pcontour(hCapprox,pwgts,zeta0,maskf,fg_run0,\
                    plot_fname=plot_fname)
    
        pcError = pcA - pcAapprox
        print('----')
        print('mean p-contour error (zeta=' + str(zeta0) + '): ' \
                + str(np.mean(np.abs(pcError))))
        print('max p-contour error (zeta=' + str(zeta0) + '): ' \
                + str(np.max(np.abs(pcError))))
        print('----')
    
    
    # hAapprox : approximation
    # Af : true
    
    # save plot 
    if save_plot:

        if map_type=='zeta':

            # plot error on zeta-contour
            plot_title = 'zeta contour error, p = ' + str(p0)
            plot_fname = os.path.join(savedir,\
                                      'zetacontour_error_' + str(p0) + '.png')
            ncontours = 21
            v = np.abs(zcError).max()
            contour_levels = np.linspace(-v,v,ncontours)

            xlim0 = fg_run0.X.min() - 0.0005
            xlim1 = fg_run0.X.max() + 0.0005
            ylim0 = fg_run0.Y.min() - 0.0005
            ylim1 = fg_run0.Y.max() + 0.0005

            formatter = FuncFormatter(long_lat)

            pl.figure();
            pl.contourf(fg_run0.X, fg_run0.Y,zcError,\
                            contour_levels,cmap='RdBu');

            # format axis
            ax = pl.gca()
            ax.yaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_formatter(formatter)
            pl.xlim([xlim0, xlim1])
            pl.ylim([ylim0, ylim1])
            for label in ax.get_xticklabels():
                label.set_rotation(30)
            pl.colorbar();
            pl.title(plot_title, size=16);
            pl.savefig(plot_fname)

            pl.close()

        elif map_type=='p':
            # plot error on p-contour
            plot_title = 'p-contour error, zeta = ' + str(zeta0)
            plot_fname = os.path.join(savedir,\
                                      'pcontour_error_' + str(zeta0) + '.png')
            ncontours = 21
            v = np.abs(np.abs(pcError).max())
            if v <= 1e-5:
                contour_levels = np.linspace(-1e-2,1e-2,ncontours)
            else:
                contour_levels = np.linspace(-v,v,ncontours)

            xlim0 = fg_run0.X.min() - 0.0005
            xlim1 = fg_run0.X.max() + 0.0005
            ylim0 = fg_run0.Y.min() - 0.0005
            ylim1 = fg_run0.Y.max() + 0.0005

            formatter = FuncFormatter(long_lat)

            pl.figure();
            pl.contourf(fg_run0.X, fg_run0.Y,pcError,\
                            contour_levels,cmap='RdBu');

            # format axis
            ax = pl.gca()
            ax.yaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_formatter(formatter)
            pl.xlim([xlim0, xlim1])
            pl.ylim([ylim0, ylim1])
            for label in ax.get_xticklabels():
                label.set_rotation(30)
            pl.colorbar();
            pl.title(plot_title, size=16);
            pl.savefig(plot_fname)

            pl.close()


    if map_type=='zeta':
        return zcError
    elif map_type=='p':
        return pcError
    else:
        print('choose map_type as "zeta" or "p"')


def save_pdata(hCapprox, pwgts, save_fname='hazardcurves_probs.npy'):
    """
    save p-contour/zeta-contour data

    """
    print('saving exceed_prob...')

    zeta1 = np.linspace(0.,  12.,1201).tolist()
    zeta2 = np.linspace(12.2,24.,  60).tolist()
    zeta = np.array(zeta1 + zeta2)

    nzeta = len(zeta)
    nevents = len(pwgts)

    exceed_prob = np.zeros((hCapprox.shape[0], len(zeta)))

    for k in range(nzeta):
        for j in range(nevents):
            exceed_prob[:,k] = np.where(hCapprox[:,j] > zeta[k], \
                                        exceed_prob[:,k] + pwgts[j], \
                                        exceed_prob[:,k])

    np.save(save_fname, exceed_prob)
    print('= done.')

def save_output(pc_mean_error_list,\
                pc_max_error_list,\
                zc_mean_error_list,\
                zc_max_error_list,\
                p0_list,zeta0_list,\
                output_fname='output_error.txt', savedir='_output'):
    """
    save error table 

    """
    output_fname = os.path.join(savedir,output_fname)
    print('saving output to: ' + output_fname)
    with open(output_fname, mode='w') as out_file:
        out_file.write('\np-contour error\n')
        out_file.write('zeta\t\tmean \t\tmax\n')
        for k in range(len(zeta0_list)):
            out_file.write(str(zeta0_list[k]) + \
                    '\t\t' + str(pc_mean_error_list[k]) + \
                    '\t\t' + str(pc_max_error_list[k]) + '\n')
        out_file.write('\nzeta-contour error\n')
        out_file.write('p\t\tmean \t\tmax\n')
        for k in range(len(p0_list)):
            out_file.write(str(p0_list[k]) + \
                    '\t\t' + str(zc_mean_error_list[k]) + \
                    '\t\t' + str(zc_max_error_list[k]) + '\n')


def read_output_Mw(Mw_list,data_dict):
    for Mw in Mw_list:
         data_dict[Mw] = read_outputs(Mw)

def reconstruct_Mw(Mw_list,data_dict,hCapprox_list,C_list, M=4):
    data_dict['Mw_list'] = Mw_list
    data_dict['fine_list'] = []
    for Mw in Mw_list:
        print('running reconstruction for magnitude = ' + str(Mw))
        hCapprox, C = run_reconstruction(Mw,data_dict,M=M)

        hCapprox_list.append(hCapprox)
        C_list.append(C)


if __name__=="__main__":

    driver_home = '../geoclaw_driver'
    savedir = '_output'
    os.system('mkdir -p ' + savedir)
    
    prb_wgts_ffname = os.path.join(driver_home,'scenario_prb_wgts.txt')
    pwgts = np.loadtxt(prb_wgts_ffname)

    # make this more flexible 
    pwgts = np.concatenate([pwgts*0.3, pwgts*0.3, pwgts*0.3, pwgts*0.1])
    
    Mw_list = [8.6, 8.8, 9.0, 9.2]
    
    # read all outputs, obtain masks, perform SVD
    print('reading all outputs')
    pkl_fname = os.path.join(savedir, 'data.pkl')
    if os.path.exists(pkl_fname):
        print('loading pickle...')
        import pickle
        with open(pkl_fname, 'rb') as handle:
            data_dict = pickle.load(handle)
    else:
        import pickle
        data_dict = {}
        read_output_Mw(Mw_list,data_dict)
        with open(pkl_fname, 'wb') as handle:
            pickle.dump(data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    # save correlation plots
    if 1:
        plot_fname = os.path.join(savedir,'cvf_Usim.png')
        fg_run0,Uc,sc,Vc,Ac,maskc,Bc,Uf,sf,Vf,Af,maskf,Bf = data_dict[9.0]
        compute_corr(Uf,Uc,M=6,plot_fname=plot_fname,save_plot=True)
           
        Yc = np.dot(np.diag(sc[:6]),Vc[:6,:]) 
        Yc = Yc / np.linalg.norm(Yc,axis=1,keepdims=True)
        Yf = np.dot(np.diag(sf[:6]),Vf[:6,:])
        Yf = Yf / np.linalg.norm(Yf,axis=1,keepdims=True)
        
        plot_fname = os.path.join(savedir,'cvf_Ysim.png')
        compute_corr(Yc.T,Yf.T,M=6,plot_fname=plot_fname,save_plot=True)

    # build approximation
    print('building approximation')
    hCapprox_list = []
    C_list = []
    reconstruct_Mw(Mw_list,data_dict,hCapprox_list,C_list,M=4)

    for k in range(len(C_list)):
        fine_list0 = data_dict['fine_list'][k]
        hCapprox_list[k][:,fine_list0] = C_list[k][:,fine_list0]


    # collect results across magnitudes
    hCapprox = np.concatenate(hCapprox_list,axis=1)
    C = np.concatenate(C_list,axis=1)
    
    S = (88,78)     # fgmax-grid shape
    maskf = np.ones(S,dtype=bool)

    fg_run0 = data_dict[Mw_list[-1]][0]
    
    zeta0_list = [0., 1., 2., 3., 4., 6., 8., 10., 12.]
    p0_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    
    pc_mean_error_list = []
    pc_max_error_list = []
    zc_mean_error_list = []
    zc_max_error_list = []
    
    # generate hazard maps
    for zeta0 in zeta0_list:
        p0 = 0.1
        pcError =\
            create_hazard_maps(hCapprox,C,pwgts,fg_run0,maskf,\
                        zeta0 = zeta0, p0=p0, map_type='p',savedir=savedir)
    
        pc_mean_error_list.append(np.mean(np.abs(pcError))) 
        pc_max_error_list.append(np.max(np.abs(pcError))) 
    
    for p0 in p0_list:
        zeta0 = 1.
        zcError =\
            create_hazard_maps(hCapprox,C,pwgts,fg_run0,maskf,\
                      zeta0 = zeta0, p0=p0, map_type='zeta', savedir=savedir)
    
        zc_mean_error_list.append(np.mean(np.abs(zcError))) 
        zc_max_error_list.append(np.max(np.abs(zcError))) 

    save_output(pc_mean_error_list,\
                pc_max_error_list,\
                zc_mean_error_list,\
                zc_max_error_list,\
                p0_list,zeta0_list,savedir=savedir)

    save_fname = os.path.join(savedir,'hazardcurves_probs.npy')
    save_pdata(hCapprox, pwgts, save_fname=save_fname)


