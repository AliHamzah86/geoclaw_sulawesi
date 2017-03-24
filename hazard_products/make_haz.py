import os,sys
from numpy import integer, size
from numpy import array, zeros, linspace, shape, nan, nanmax
from numpy import ma, where
from numpy import load, save, loadtxt

## Not used at the moment
def combine_prob(p1,p2):
    """Returns the probability that event 1 or 2 happens"""
    return 1. - (1-p1)*(1-p2)


def hazard_curves_calcv(hazard_curves_file,fixed_grid_file,events_directories,events_probs,\
                        zeta,fgno,comments_file=None):
    
    # Note: since the fixed_grid_file will be x,y,B0, there will only be one for all
    #       the events for this fgno, but its dimensions will change depending on the fgno.

    d=load(fixed_grid_file)
    sh=shape(d)
    npoints=sh[0]

    # Store the initial Bathymetry in B0. This is relative to MHW.
    xlong=d[:,0]; ylat=d[:,1]; B0 = d[:,2];               

    #masking = where(B0 < 0,True,False)

    print 'INSIDE ROUTINE hazard_curves_calcv:'
    print ' '
    print 'The fixed grid number was: ',fgno
    print 'The number of points on the fixed grid was: ',npoints 
    print ' '
    print 'The zeta values inside make_haz were: '
    print zeta
    print ' '

    #Initialize the hazard curves file called exceed_prob and mask out
    #where the initial bathy B0 is negative. No, now not masking, computing zeta=h+B0
    #for points where B0 is negative.
    #
    nzeta = len(zeta)
    exceed_prob = zeros((npoints,nzeta))                    # initialize to zero

    #exceed_prob = ma.masked_array(zeros((npoints,nzeta)))  # initialize to zero
    #for iz in range(nzeta):
    #    exceed_prob[:,iz]=ma.masked_where(masking,exceed_prob[:,iz])

    # loop over all events and update exceed_prob at each grid point by combining
    # current value with the probability Pk of this event:
    #
    nevents=len(events_directories)
    for event in range(nevents):
        event_dir=events_directories[event]
        #event_dir = os.path.join(events_dir, event)
        fgname = 'hmax_fg' + str(fgno) + '.npy'
        hmax_file = os.path.join(event_dir, fgname)
        print 'PROCESSING the hmax_fg file: '
        print hmax_file
        print ' '
        hmax = load(hmax_file)                    #contains h of Geoclaw, not zeta        
                                                  #so for deeper water, hmax really big
                                                  #so work with zeta = hmax + B0 there
        #### debugging
        #print 'In make_haz, For event: ',event
        #print 'hmax[0:100] after initial loading'
        #print hmax[0:100]
        #### debugging

        waterpt=B0 < 0.0
        hmax = where(waterpt, hmax+B0,hmax)       #hmax now is h on land, h+B0 where B<0
                                                  #we used to call this zeta
                                                  #hmax could potentially be negative with roundoff
        #### debugging
        #print ' '
        #print 'hmax[0:100] after adding B0 for waterpoints'
        #print hmax[0:100]
        #print ' '
        #### debugging

        #hmax = ma.masked_where(masking,hmax)      #Question: Do I have to mask hmax?
        #Hmax = hmax.reshape((nx,ny),order='F')

        for k in range(nzeta):
            Pk = exceed_prob[:,k]  # conditional probabilities at all points for one 
                                   # exceedance value zeta_k
            #exceed_prob[:,k] = where(hmax > zeta[k], \
            #                   combine_prob(events_probs[event],Pk), Pk)
            exceed_prob[:,k] = where(hmax > zeta[k], \
                               events_probs[event]+Pk, Pk)
        
    ##### debugging
    print 'The rows 0:100 of the first column of exceed_prob (for zeta=0) was: '
    print exceed_prob[0:100,0]
    print ' '
    ##### debugging

    print ' ' 
    print 'Computed exceedance probabilities for fgno= ',fgno
    print ' '
    print 'The probabilities are stored in file: '
    print hazard_curves_file
    print ' '

    # .npy files can't have mask, so replace with nans where it is masked before
    # storing as an .npy file
    #
    #for iz in range(nzeta):
    #    exceed_prob[:,iz]=where(masking,nan,exceed_prob.data[:,iz])

    print 'Maximum exceedance over all grid points is %8.5f' % nanmax(exceed_prob)
    print ' '

    save(hazard_curves_file,array(exceed_prob))

    print 'FIRST GRID POINT OUTPUT:'
    print 'The first grid point in exceed_prob had long, lat, B0: '
    print xlong[0], ylat[0], B0[0]
    print ' '
    print 'The nzeta exceed_prob values after mask=nan were: '
    for iz in range(nzeta):
        print zeta[iz],exceed_prob[0,iz]
    print ' '
    print 'LAST GRID POINT OUTPUT:'
    print 'The last grid point in exceed_prob had long, lat, B0: '
    print xlong[-1], ylat[-1], B0[-1]
    print ' '
    print 'The nzeta exceed_prob values after mask=nan were: '
    for iz in range(nzeta):
        print zeta[iz],exceed_prob[-1,iz]
    print ' '
    #


if __name__ == '__main__':

    projectdir = os.environ['FEMA']
    scenarioname = 'scenario_test2'
    scenariodir = projectdir + '/' + scenarioname
    scenario_weights_file=scenariodir + '/scenario_prb_wgts.txt'
    runs_todo_file = scenariodir + '/runs_todo.txt'
    fg_types_todo_file = scenariodir + '/fg_types_todo.txt'
    fg_transects_todo_file = scenariodir + '/fg_transects_todo.txt'
    fg_2Dgrids_todo_file = scenariodir + '/fg_2Dgrids_todo.txt'

    print 'INSIDE if name is main, data echo: '
    print ' '
    print 'The project directory was: ',projectdir
    print ' '
    print 'The scenario directory was: ',scenariodir
    print ' '
    
    runs_todo = loadtxt(runs_todo_file,dtype=integer)
    print 'The master run numbers used for hazard curves are: '
    print runs_todo
    print ' '

    fg_types_todo = loadtxt(fg_types_todo_file,dtype=integer)
    print 'fg_types_todo where 1 means transects, 2 means 2Dgrids is: '
    print fg_types_todo
    print ' '

    events_directories=[]
    for event in runs_todo:
        events_directories.append(projectdir+'/fine_runs/run_'+str(event))
    print 'The directories for all the events was: '
    print events_directories 
    print ' '

    events_probs = loadtxt(scenariodir + '/scenario_prb_wgts.txt')
    print 'The conditional probability weights for scenario events was: '
    print events_probs
    print ' '

    zeta = linspace(0,12,1201)
    print 'The 1201 exceedance values were zeta=linspace(0,12,1201)'
    print ' '

    ######  Now loop over all the work to be done, and call hazard_curves_calcv ###
    ######  for each fixed grid of each type.
    for itype in fg_types_todo:
        if (itype == 1):            #Do all the transect type of fixed grids
            fg_transects_todo=loadtxt(fg_transects_todo_file,dtype=integer)
            if (size(fg_transects_todo) == 1):
                fg_transects_todo = array([fg_transects_todo])
            for iwhich in fg_transects_todo:
                outdir = scenariodir + '/' + scenarioname + '_transect_fg' + str(iwhich)
                hazard_curves_file = outdir + '/hazardcurves_probs.npy'
                fixed_grid_file = projectdir+'/fine_runs/'+'fine_xyB0_fg'+str(iwhich) +'.npy'
                print '#####################'
                print 'CALLING hazard_curves_calcv with fixed grid number: ',iwhich
                print 'The fixed grid xyB0 file was: ',fixed_grid_file
                print ' '
                print 'The hazard curves .npy file was: '
                print hazard_curves_file
                print ' '
                hazard_curves_calcv(hazard_curves_file,fixed_grid_file,\
                                    events_directories,events_probs,\
                                    zeta,iwhich,comments_file=None)
                print ' '
        if (itype == 2):            #Do all the 2D grid type of fixed grids
            fg_2Dgrids_todo=loadtxt(fg_2Dgrids_todo_file,dtype=integer)
            if (size(fg_2Dgrids_todo) == 1):
                fg_2Dgrids_todo = array([fg_2Dgrids_todo])
            for iwhich in fg_2Dgrids_todo:
                outdir = scenariodir + '/' + scenarioname + '_2Dgrid_fg' + str(iwhich)
                hazard_curves_file = outdir + '/hazardcurves_probs.npy'
                fixed_grid_file = projectdir+'/fine_runs/'+'fine_xyB0_fg'+str(iwhich) +'.npy'
                print '#####################'
                print 'CALLING hazard_curves_calcv with fixed grid number: ',iwhich
                print 'The fixed grid xyB0 file was: ',fixed_grid_file
                print ' '
                print 'The hazard curves .npy file was: '
                print hazard_curves_file
                print ' '
                hazard_curves_calcv(hazard_curves_file,fixed_grid_file,\
                                    events_directories,events_probs,\
                                    zeta,iwhich,comments_file=None)
                print ' '
