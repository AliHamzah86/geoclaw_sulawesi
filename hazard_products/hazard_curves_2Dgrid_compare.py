"""

Usage: 
 >>>
 from hazard_curves_2Dgrid_compare import hazard_curves_2Dgrid_compare

    where for example

    INPUT: 
      hazard_curves_file = the .npy file used for the plots:

           P1 P2 P3 PNzeta 

      fixed_grid_file = file with x y B in colms, a .npy file
        (for the FEMA project, B is B0, initial bathy)
        Made a list, so the scenarios being compared could each
        have a different B value, say to compare a fine with a 
        coarse scenario.  The x and y are assumed the same. To
        tell which scenario is which, we added the list
        scenarioresolution below

      scenarioresolution = ['fine','fine'] to compare two fine ones,
                         = ['fine','coarse'] when first is fine, 2nd coarse.

      zeta = numpy array of Nzeta exceedance values 

      whichrows = numpy array of integers representing the
                  rows of the hazard_curves_file to plot
                  the first row in numpy array is index 0!
                  Each row represents one hazard curve at a
                  particular point.

      whichzetas = numpy array of integers representing the
                   columns in zeta for which we want contour
                   plots of probabilities. For example, 
                   whichzetas=array([0,3]) would use the
                   probabilities corresponding to P1(col 0),
                   and P4(col 3).  Each of the whichzetas is
                   an exceedance value (of depth). So
                   plotting what is in a column is a p-contour.
                   For example a representation
                   of the probabilities at all the points of
                   exceeding a whichzeta value (say 2 meters).
                   

      Eventpath  = directory path where to save the plots

      pbarlist: 
             Made this input, rather than set within the code.  It is
             a list of "floods" for the zeta-contours.
             For example, pbarlist=[.01,.002] would be the 100 year flood
             or the 500 year flood.  Another example, if we are only interested
             in conditional probabilities (not whether or not Cascadia happens,
             for example), we might take pbarlist=[.3,.5,.7,.8,.9].  The
             entries are the probability that a level is exceeded. 

      Npts1, Npts2:  The 2D fixed grid has Npts1 points in the longitude and Npts2
                     points in the latitude direction.  These are used to reshape
                     data for plotting.

      X_CC,Y_CC,Z_CC:
             These are 2D structures for the little contour plot in upper right
             corner of some plots.

      legendnames:  names of the scenarios used in the comparison to be put in a
                    legend on the hazard curves Comparison plots.


      The following region_extent will be calculated in the routine.

      region_extent: [longwest,longeast,latsouth,latnorth] for the 2D fixed grid. Use
                     for plotting.


     Example:
         import os
         CCdir = os.environ['CCdir']
         probsnpy='hazardcurves_probs.npy'
         Testsdir1= CCdir + '/Realizations_Test1/'
         Testsdir2= CCdir + '/Realizations_Test2/'
         fixed_grid_file = CCdir + '/Datafiles/fixedgrid_xyB_test.npy'

         hazard_curves_file  = [Testsdir1 +'EVENTeall' + '/'+ probsnpy,
                               [Testsdir2 +'EVENTeall' + '/'+ probsnpy]

         OutZetafile = [Testsdir1 + 'EVENTeall' + '/' + 'Zeta100_500.npy',
                        Testsdir2 + 'EVENTeall' + '/' + 'Zeta100_500.npy]

         Eventpath = CCdir + 'Realizations/comparisons/comparisons_1/2Dgrid_fg1'

         See the __main__ below for setting zeta

         speed=True if we are doing plots of speeds; otherwise speed=False.
         depth=True if we are doing plots of depth; otherwise depth=False.
         mflux=True if we are doing plots of mflux; otherwise mflux=False.

    OUTPUT:  Hazard Curves comparison plots for the whichrows grid points.
             p-contour comparisons plots for various zeta levels in whichzetas.
             zeta-contour comparisons plots for various p levels in pbarlist.

             These comparison plots are all stored in the directory Eventpath
             (which is called comparisondir in the calling program). A commentsfile
             of the job run is also stored there.



"""
from numpy import loadtxt
from numpy import array,shape,load,nan,isnan,linspace,sum,where
from numpy import zeros,ones,ma,amin,amax,logical_or,logical_and
from numpy import int,integer,size,floor,ceil
from numpy import cos,pi,hstack,arange
from numpy import nanmax,nanmin,nanmean
from matplotlib.mlab import find
from matplotlib import image,colors
import matplotlib.pyplot as plt
import os,sys

def hazard_curves_2Dgrid_compare(hazard_curves_file,fixed_grid_file,\
       OutZetafile,\
       zeta,whichzetas,Eventpath,pbarlist,Npts1,Npts2,\
       X_CC,Y_CC,Z_CC,legendnames,scenarioresolution,\
       whichrows=None,speed=False,depth=True,mflux=False):

    #Specific to Crescent City and the FEMA project
    #FEMA=os.environ['FEMA']                           #If we need this
    #CCmap=image.imread(FEMA+'/maps/CCmap8.png')       #don't have this yet

    aspectNO=41.75                                     #aspect latitude for plotting
                                                       #specific to Crescent City

    commentsfile=Eventpath + '/hc_compare_run.output'
    outfile=open(commentsfile,'w')

    outfile.write(' THE 2D fixed grid dimensions Npts1 and Npts2 are = : ')
    outfile.write('%s %s ' %(Npts1,Npts2))
    outfile.write(' \n')
    outfile.write(' \n')

    #Load the hazard curves .npy files for each scenario and put those arrays in a list
    num_scenarios=len(hazard_curves_file)
    scenario_rows=int(ceil(num_scenarios/2.0))   #for organizing subplots
    print 'NUM_SCENARIOS was: ',num_scenarios
    print 'SCENARIO_ROWS was: ',scenario_rows

    plt.figure(10,(12,8))                      #compare the bathymetry files
    fprobs_list=[]; zetapbar_list=[]; B_list=[]; B0_pos_list=[];
    for j in range(num_scenarios):
        outfile.write(' \n')
        outfile.write(' Scenario: %s ' %legendnames[j])
        outfile.write(' \n')
        d=load(fixed_grid_file[j])
        x=d[:,0]; y=d[:,1];                       #long and lat for this scenario
        B=d[:,2]                                  #This is B0 for this scenario for the FEMA project
        B_list.append(B)
        B0_pos_list.append(where(B < 0.0, 0, B))  #Used to make eta. If water point, B already added
        fprobs_list.append(load(hazard_curves_file[j]))
        zetapbar_list.append(load(OutZetafile[j]))

        plt_number=j+1
        plt.subplot(scenario_rows,2,plt_number)
        Bsquare=B.reshape((Npts1,Npts2),order='F') 
        xsquare=x.reshape((Npts1,Npts2),order='F')
        ysquare=y.reshape((Npts1,Npts2),order='F')
        plt.contour(xsquare,ysquare,Bsquare,arange(0,21,2), colors='k')
        plt.ticklabel_format(format='plain',useOffset=False)
        plt.xlim(x[0],x[-1])
        if (logical_or(plt_number == num_scenarios,plt_number == num_scenarios-1)):
            plt.xticks(rotation=20,fontsize=10)
            plt.yticks(fontsize=10)
        else:
            plt.xticks([])
            plt.yticks(fontsize=10)
        plt.ylim(y[0],y[-1])
        a = plt.gca()
        a.set_aspect(1./cos(aspectNO*pi/180.))
        titlestr= '\n' + 'B-Contours for  ' + legendnames[j]
        plt.title(titlestr)
    outfile.write(' \n')

    fname=Eventpath + '/bathymetry-contours.png'
    plt.savefig(fname)
    plt.close(10)

    #All the scenarios have the same x and y, so use the last one from above that is in x and y

    #Compute the region_extent from the x and y information for this 2Dgrid fixed grid
    region_extent=[]; region_extent.append(x[0]); region_extent.append(x[-1]);
    region_extent.append(y[0]); region_extent.append(y[-1]);
    outfile.write(' THE region_extent for this fgno fixed grid is : \n')
    outfile.write('%s ' %region_extent)
    outfile.write(' \n')
    outfile.write(' \n')

    if (depth==True):
        outfile.write(' This is a DEPTH zeta run.  \n')
        outfile.write(' \n')
        outfile.write(' Comparison plots in directory: \n')
        outfile.write('  %s ' %Eventpath)
        outfile.write(' \n')
        outfile.write(' \n')
        outfile.write(' Setting mflux and speed to False')
        outfile.write(' \n')
        mflux=False
        speed=False

    else:
        print ('program stopping, depth was not True')
        sys.exit(1)

    outfile.write('  \n')
    outfile.write(' The zeta_i exceedance values follow. \n')
    outfile.write('%s ' %zeta)
    outfile.write(' \n')
    outfile.write(' \n')

#   Compute the rows (points) for which a hazard curve plot will be made.
#   We have Npts1 grid points in the long-direction and Npts2 in the lat-direction.
#   Cover this grid with 5 hazard curves in each direction, for a total of 25.
    
#   Say we want the first row of grid points to be taken from row jj=int(Npts2/5.0).
#   Then we skip the (jj-1)*Npts1 grid points in the (jj-1) rows below.  Then the
#   points in this row jj will start in column ii=int(Npts1/5.0) - 1, and continue to
#   column 5*ii -1, the -1 due to python starting at 0.

    if (whichrows == None):
        whichrows=[]

    ii=int(Npts1/5.0); jj=int(Npts2/5.0);
    for rows in [jj,2*jj,3*jj,4*jj,5*jj]:
        for cols in [ii,2*ii,3*ii,4*ii,5*ii]:
            whichindex=(rows-1)*Npts1 + cols -1    #python arrays begin with 0
            whichrows.append(whichindex)
    whichrows=array(whichrows)
    outfile.write(' The location indices for hazard plots are: \n')
    outfile.write('%s ' %whichrows)
    outfile.write(' \n')
    outfile.write(' \n')

    outfile.write(' The indices for zeta-exceedance levels for p-contours are: \n')
    outfile.write('%s ' %whichzetas)
    outfile.write(' \n')
    outfile.write(' \n')

    outfile.write(' The flood levels (p values) for zeta-contours are: \n')
    outfile.write('%s ' %pbarlist)
    outfile.write(' \n')
    outfile.write(' \n')

    #getting the x and y ready for pcolor or contourf plots as xp and yp
    #Note: We only need to do this is the fixed grid is of type fgno=3;
    #
    #B=B.reshape((Npts1,Npts2),order='F')    #choosing from B_pos_list and shaping later
    #
    dx=.000001/2.; dy=dx;
    xp=zeros((Npts1+1,Npts2+1))
    yp=zeros((Npts1+1,Npts2+1))
    xp[0:Npts1,0:Npts2]=x.reshape((Npts1,Npts2),order='F')
    xp[0:Npts1,Npts2]=xp[0:Npts1,Npts2-1]        #add another colm for same x-value
    xp[Npts1,:]=xp[Npts1-1,:]+dx                 #add another row with bigger x-value
    yp[0:Npts1,0:Npts2]=y.reshape((Npts1,Npts2),order='F')
    yp[0:Npts1,Npts2]=yp[0:Npts1,Npts2-1]+dy     #add another col with bigger y-value
    yp[Npts1,:]=yp[Npts1-1,:]                    #add another row for same y-value           

    ### Now do the hazard curves Comparison Plots
    plt.figure(11,(12,8))
    fprobs121=zeros(121)        
    for irow in whichrows:       #Grid Point Hazard Curve Plots
        plt.clf()
        B_vals=[]                #Accumulate for the plot title, a check too
        for j in range(num_scenarios):
            B_vals.append(B_list[j][irow])
            fprobs=fprobs_list[j]
            fprobs121=fprobs[irow,:]
            masking=where(isnan(fprobs121),True,False)
            fprobs121=ma.masked_where(masking,fprobs121)

            #fprobs121=where(fprobs121 <= 1e-10,1e-10,fprobs121) #So semilog can work

            if (mflux == True):
                fprobs121=where(fprobs121 <= 1e-10,1e-10,fprobs121) #mflux use loglog
                plt.loglog(zeta,fprobs121,'o-',markersize=5)

            else:                                    #Not doing semilog for depth or speed plots 
                                                     #We don't have recurrence time of earthquake
                                                     #We just have conditional probs of realizations
                #plt.plot(zeta,fprobs121,'o-',markersize=5)
                plt.plot(zeta,fprobs121,'-')

                #plt.semilogy(zeta,fprobs121,'ro')
                #plt.hold(True) 
                #plt.semilogy(zeta,fprobs121,'r-')
            plt.hold(True)
        #ylim(0,.040)
        #ylim(0,.010)
        plt.xlim(0.0,zeta[-1])
        plt.ylim(-0.1,1.1)                  
        plt.legend(legendnames,loc='center right')

        #Not using semilogy so comment this out and choose ylim automatically
        #plt.ylim(1e-10,100)

        #Use yticks for semilogy plots
        #yticks([1e-5,1e-4,1e-3,.002,.01,.1,1],\
        #       ['1e-5','1e-4','1e-3','0.002','0.01','0.1','1'])
        #yticks([1e-10,1e-7,1e-6,1e-5,1e-4,.0004,1e-3,.002,.01,.1,1],\
        #  ['1e-10','1e-7','1e-6','1e-5','1e-4','.0004','1e-3','0.002','0.01','0.1','1'])

        #Use yticks for tideprobs, or for the conditional probs of realizations
        plt.xticks([0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0],\
           ['0','2','4','6','8','10','12','14','16','18','20','22','24'])
        plt.yticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],\
                   ['0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'])

        nn=str(irow).zfill(6)
        if (depth==True):
            plt.xlabel('Exceedance Level for Zeta (m)')
            fname=Eventpath + '/HazardCurve_Zeta_GridPoint_%s.png' %nn
        elif (speed==True):
            plt.xlabel('Exceedance Level for Speed (m/sec)')
            fname=Eventpath + '/HazardCurve_Speed_GridPoint_%s.png' %nn
        else:
            plt.xlabel('Exceedance Level for Momentum Flux ')
            fname=Eventpath + '/HazardCurve_Mflux_GridPoint_%s.png' %nn

        plt.ylabel('Probability of Exceedance')
        plt.title('At Longitude, Latitude (%9.4f,%9.4f) \n' % (x[irow],y[irow]) \
                  +'Scenario Bathymetry Values: ' + num_scenarios*'%9.4f '\
                   %tuple(B_vals))


        ######### putting little map on the plot, showing the location
        #fixed grid was a 2Dgrid
        plt.axes([.7,.65,.25,.25])  #.7,.65 is lower left ratio of figure,.25,.25 wdth,height
        plt.contour(X_CC, Y_CC, Z_CC, arange(0,21,2), colors='k')
        a = plt.gca()
        a.set_aspect(1./cos(aspectNO*pi/180.))
        xi=x[irow]; yi=y[irow];
        plt.plot([xi],[yi],'ro')
        plt.xticks([])
        plt.yticks([])
        #########

        plt.savefig(fname)    #save plot for this irow in whichrow
        #raw_input(' hit return to get the next plot ')
    plt.close(11)

    ### Now do the p-contour Comparison Plots
    plt.figure(15,(12,8))         #one of these for each whichzeta
    for izeta in whichzetas:      #p-contour plots for heights in whichzetas
        plt.clf()                 #clears figure 15
        for j in range(num_scenarios):
            fprobs=fprobs_list[j]
            zetarest=fprobs[:,izeta].reshape((Npts1,Npts2),order='F')
            masking=logical_or(isnan(zetarest),zetarest <=0.0)
            # masking is where zetarest is masked or a probability was 0 that 
            # level was exceeded, so will remain white on the graph.
            # masking=isnan(zetarest)
            zetarest=ma.masked_where(masking,zetarest)

            #raw_input(' hit return for p-contour for zeta= %5.2f ' %(zeta[izeta]))
            #As the j-loop progresses, plot the contours for the scenarios left to
            #right top to bottom, using 2 columns in each row of plots
            plt_number=j+1
            plt.subplot(scenario_rows,2,plt_number) 

            #raw_input(' hit return for p-contour for zeta= %5.2f ' %(zeta[izeta]))
            #V=[0.0,.0004,.001,.0015,.002,.0025,.005,.0075,.01,.015,.02,1.]
            #V=[0.0,.0004,.001,.002,.004,.005,.01,.02,.04,.06,1.]
            # The V below is for the Discovery Bay, Lynch Cove tideprob runs
            # or for the FEMA project which has conditional probabilities
            V2=[0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.000001]
            V=[0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
            #0,2500yr,1000yr,500yr,250yr,200yr,100yr,50yr,25yr,20yr,always
            #V=[0.0,.00001,.0001,.0002,.0004,.0008,.0010,.0012,.0014,.0016,.0018,.0020,.0040]
            norm=colors.BoundaryNorm(V,256)
            CS=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],zetarest,V2,norm=norm)
            if (num_scenarios == 2):
                cbar_shrink=.7
            else:
                cbar_shrink=.9
            cbar = plt.colorbar(shrink=cbar_shrink)
            cbar.set_ticks(V)
            plt.hold(True)

            # Contours of topo added:
            if 0:   #Would need to choose the right B and reshape for this one
                clines2 = linspace(0,20,11)
                plt.contour(xp[:-1,:-1],yp[:-1,:-1],B,clines2,colors='k')
                plt.ticklabel_format(format='plain',useOffset=False)
                plt.xticks(rotation=20)
                plt.axis(region_extent)
                a = plt.gca()
                a.set_aspect(1./cos(aspectNO*pi/180.))
            if 1:
                #plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
                plt.contour(X_CC,Y_CC,Z_CC, arange(0,8,4), colors='k')
                #Note: only plt xlim if plt_number == num_scenarios or if
                #      plt_number == num_scenarios-1
                #      Then the xlim only occurs in the last row of plots
                #plt.xlim(-124.215,-124.18)
                plt.xlim(region_extent[0],region_extent[1])
                if (logical_or(plt_number == num_scenarios,plt_number == num_scenarios-1)):
                    plt.xticks(rotation=20,fontsize=10)
                    plt.yticks(fontsize=10)
                else:
                    plt.xticks([])
                    plt.yticks(fontsize=10)
                #plt.ylim(41.735,41.768)
                plt.ylim(region_extent[2],region_extent[3])
                plt.ticklabel_format(format='plain',useOffset=False)

            if 0:
                #Draw green box around the xp,yp contour sub-plot
                box_x=[region_extent[0],region_extent[1],region_extent[1],\
                       region_extent[0],region_extent[0]]
                box_x=array(box_x)
                box_y=[region_extent[2],region_extent[2],region_extent[3],\
                       region_extent[3],region_extent[2]]
                box_y=array(box_y)
                plt.plot(box_x,box_y,'g-',linewidth='2')

            a = plt.gca()
            a.set_aspect(1./cos(aspectNO*pi/180.))

            nn=str(int(zeta[izeta]*100)).zfill(5)
            if (depth==True):
                plt.title('\np-contours, zeta=%5.2f m, %s ' %(zeta[izeta],legendnames[j]), fontsize=10)
                fname=Eventpath + '/p-contours_zeta_%s.png' %nn
            elif (speed==True):
                plt.title('p-contours Comparisons for speed=%5.2f m/sec' %(zeta[izeta]), fontsize=10)
                fname=Eventpath + '/p-contours_speed_%s.png' %nn
            else:
                plt.title('p-contours Comparisons for momentum flux=%5.2f ' %(zeta[izeta]), fontsize=10)
                fname=Eventpath + '/p-contours_mflux_%s.png' %nn
        plt.savefig(fname)
        plt.clf()

        outfile.write('\n')
        outfile.write(' Processing p-contour diffs for zeta= %5.2f \n' \
                   %zeta[izeta])
        fprobs0=fprobs_list[0]
        zetarest0=fprobs0[:,izeta].reshape((Npts1,Npts2),order='F')    #not masked yet
        for j in range(num_scenarios):
            fprobs=fprobs_list[j]
            zetarest=fprobs[:,izeta].reshape((Npts1,Npts2),order='F') #not masked yet

            diffplot=abs(zetarest-zetarest0)
            masking=logical_or(isnan(diffplot),diffplot <=0.0)
            # masking is where probability differences are 0.0. Dont think
            # diffplot will have any nans in it at all, since fprobs is not masked
            # for water points now and the probs will be 0 to 1. 
            # If the prob differences are 0, want to remain white on the graph.
            probdiffmax=nanmax(diffplot)
            #diffplot_pos=where(diffplot > 0.0,diffplot,nan)  ##not for the probabilities
            #probdiffmin=nanmin(diffplot)
            probdiffmean=nanmean(diffplot)
            outfile.write('\n')
            outfile.write('max value in probdiffplot in %s was  %9.5f \n' \
                      %(legendnames[j],probdiffmax))
            #outfile.write('min value in probdiffplot in %s was  %9.5f \n' \
            #         %(legendnames[j],probdiffmin))
            outfile.write('mean of probdiffplot in %s was  %9.5f \n' \
                      %(legendnames[j],probdiffmean))
            ##Now do the masking for the plot
            diffplot=ma.masked_where(masking,diffplot)

            plt_number=j+1
            plt.subplot(scenario_rows,2,plt_number) 

            # The V below is for the Discovery Bay, Lynch Cove tideprob runs
            # or for the FEMA project which has conditional probabilities
            V2=[0.0,.025,.05,.075,.1,.2,.3,.4,.6,.8,1.000001]
            V=[0.0,.025,.05,.075,.1,.2,.3,.4,.6,.8,1.]
            norm=colors.BoundaryNorm(V,256)
            CS=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],diffplot,V2,norm=norm)
            if (num_scenarios == 2):
                cbar_shrink=.7
            else:
                cbar_shrink=.9
            cbar = plt.colorbar(shrink=cbar_shrink)
            cbar.set_ticks(V)
            plt.hold(True)

            # Contours of topo added:
            if 1:
                #plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
                plt.contour(X_CC,Y_CC,Z_CC, arange(0,8,4), colors='k')
                #Note: only plt xlim if plt_number == num_scenarios or if
                #      plt_number == num_scenarios-1
                #      Then the xlim only occurs in the last row of plots
                #plt.xlim(-124.215,-124.18)
                plt.xlim(region_extent[0],region_extent[1])
                if (logical_or(plt_number == num_scenarios,plt_number == num_scenarios-1)):
                    plt.xticks(rotation=20,fontsize=10)
                    plt.yticks(rotation=10)
                else:
                    plt.xticks([])
                    plt.yticks(rotation=10)
                #plt.ylim(41.735,41.768)
                plt.ylim(region_extent[2],region_extent[3])
                plt.ticklabel_format(format='plain',useOffset=False)

            if 0:
                #Draw green box around the xp,yp contour sub-plot
                box_x=[region_extent[0],region_extent[1],region_extent[1],\
                       region_extent[0],region_extent[0]]
                box_x=array(box_x)
                box_y=[region_extent[2],region_extent[2],region_extent[3],\
                       region_extent[3],region_extent[2]]
                box_y=array(box_y)
                plt.plot(box_x,box_y,'g-',linewidth='2')

            a = plt.gca()
            a.set_aspect(1./cos(aspectNO*pi/180.))

            nn=str(int(zeta[izeta]*100)).zfill(5)
            if (depth==True):
                plt.title('\np-diff-contours, zeta=%5.2f m, %s ' %(zeta[izeta],legendnames[j]), fontsize=10)
                fname=Eventpath + '/p-diff-contours_zeta_%s.png' %nn
        plt.savefig(fname)
    plt.close(15)

    #Now do the "Floods" Comparison contour plots of zeta(h) and eta(h+B0)
    titlelist=[]
    for j in range(len(pbarlist)):
        titlelist.append('p=' + str(pbarlist[j]))
    
    plt.figure(14,(12,8))                  #each pbar individually

    #HERE, saving a weighted average of zetadiff across pbars
    zetadiffplot_weighted=[]
    for jsc in range(num_scenarios):
        zetadiffplot_weighted.append(zeros((Npts1,Npts2)))
    sumpbarlist=sum(array(pbarlist))
    outfile.write(' \n')
    outfile.write(' sum of pbarlist was: %9.2f ' %sumpbarlist)
    outfile.write(' \n')
    #HERE

    for j in range(len(pbarlist)):
        plt.clf()
        outfile.write(' \n')
        outfile.write(' FOR PBAR == : %s ' %pbarlist[j])
        for jsc in range(num_scenarios):
            zetapbar=zetapbar_list[jsc]
            zetaplot=zetapbar[:,j].reshape((Npts1,Npts2),order='F')
            nanmax_zetaplot=nanmax(zetaplot)
            #nanmin_zetaplot=nanmin(zetaplot)

            #masking=logical_or(isnan(zetaplot),zetaplot <=0.0)
            #When zetaplot was a nan at B0 < 0 points and where no flood, this
            #would mask where B0<0 or no zeta from 0 to 12 exceeded (a nan)
            #and also where zeta<=0.0 was the max exceeded with this probability.
            #I am not sure 0.0 should be included above in the masking.

            masking=isnan(zetaplot)
            #This will mask where no zeta h or h+B0 from 0 to 12 was exceeded with
            #this probability, and had been stored as a nan -- no flood situation.
            zetaplot=ma.masked_where(masking,zetaplot)

            cmax = ceil(max(nanmax_zetaplot,5.)*100.)/100.
            #cmax = ceil(max(zetaplot.max(), 5.)*100.)/100.

            outfile.write('\n')
            outfile.write('max value in zetaplot in %s was  %9.5f \n' \
                      %(legendnames[jsc],nanmax_zetaplot))
            #         %(amax(zetaplot)))
            #outfile.write('min value in zetaplot in %s was  %9.5f \n' \
            #          %(legendnames[jsc],nanmin_zetaplot))
            #         %(amin(zetaplot)))

            if (depth==True):
                if (cmax <= 16.0):
                    clines = [1e-3] + [.5,1.0] + list(linspace(2.0,16.0,8))
                else:
                    clines = [1e-3] + [.5,1.0] + list(linspace(2.0,14.0,7)) + [cmax]
            elif (speed==True):
                dmax=max(cmax,10.0)
                clines = [1e-3] + list(linspace(1,9,9)) + [dmax]
            else: 
                dmax=max(cmax,1200.)
                #clines = [1e-3] + [25.0,50.0,100.0,150,200.0,300.0,400.0,600.0,800.0] +\
                          #[1000.0] + [dmax]

                clines=[0.01, 0.05]+list(linspace(0.1,2,20))+\
                       [4,6,8,10,20,50,100,200,400,600,1000,dmax]
            nlines = len(clines)
            n1 = int(floor((nlines-1)/2.))
            n2 = nlines - 1 - n1
            Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
            Red = hstack([linspace(0,0.8,n1), ones(n2)])
            Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
            colors2 = zip(Red,Green,Blue)

            #As jsc marches through the scenarios, the contour plots are plotted
            #using plots 2 columns wide, top to bottom, left to right
            plt_number=jsc+1
            plt.subplot(scenario_rows,2,plt_number) 
            CS100=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],zetaplot,\
                       clines,colors=colors2)
            if (num_scenarios == 2):
                cbar_shrink=.7
            else:
                cbar_shrink=.9
            cbar = plt.colorbar(shrink=cbar_shrink)
            cbar.set_ticks(clines)
            plt.hold(True)

            # Contours of topo added:
            if 0:  #Must choose and fix the B if this is to be used
                clines2 = linspace(0,20,11)
                plt.contour(xp[:-1,:-1],yp[:-1,:-1],B,clines2,colors='k')
                plt.ticklabel_format(format='plain',useOffset=False)
                plt.xticks(rotation=20)
                plt.axis(region_extent)
                a = plt.gca()
                a.set_aspect(1./cos(aspectNO*pi/180.))
            if 1:
                #plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
                plt.contour(X_CC,Y_CC,Z_CC, arange(0,8,4), colors='k')
                plt.ticklabel_format(format='plain',useOffset=False)
                #plt.xlim(-124.215,-124.18)
                plt.xlim(region_extent[0],region_extent[1])
                if (logical_or(plt_number == num_scenarios,plt_number == num_scenarios-1)):
                    plt.xticks(rotation=20,fontsize=10)
                    plt.yticks(fontsize=10)
                else:
                    plt.xticks([])
                    plt.yticks(fontsize=10)
                #plt.ylim(41.735,41.768)
                plt.ylim(region_extent[2],region_extent[3])

            if 0:
                #Draw green box around the xp,yp sub-contour plot
                box_x=[region_extent[0],region_extent[1],region_extent[1],\
                       region_extent[0],region_extent[0]]
                box_x=array(box_x)
                box_y=[region_extent[2],region_extent[2],region_extent[3],\
                      region_extent[3],region_extent[2]]
                box_y=array(box_y)
                plt.plot(box_x,box_y,'g-',linewidth='2')

            a = plt.gca()
            a.set_aspect(1./cos(aspectNO*pi/180.))

            if (depth==True):
                titlestr='\nzeta-contours, p= ' + str(pbarlist[j]) + ', ' + legendnames[jsc]
                plt.title(titlestr,fontsize=10)
                fname=Eventpath + '/zeta-contours_prob_' + titlelist[j] + '.png'
            elif (speed==True):
                titlestr='speed-contours Comparisons for p= ' + str(pbarlist[j])
                plt.title(titlestr)
                fname=Eventpath + '/speed-contours_prob_' + titlelist[j] + '.png'
            else:
                titlestr='mflux-contours Comparisons for p= ' + str(pbarlist[j])
                plt.title(titlestr)

        fname=Eventpath + '/zeta-contours_prob_' + titlelist[j] + '.png'
        plt.savefig(fname)
        plt.clf()        #Clear figure 14 after all scenario plots are on the figure and saved

        ##Finding the difference of zetapbar for all scenarios from zetapbar of scenario jsc=0 
        zetapbar0=zetapbar_list[0]
        zetaplot0=zetapbar0[:,j].reshape((Npts1,Npts2),order='F')
        nan0=isnan(zetaplot0)
        notnan0 = ~isnan(zetaplot0)
        for jsc in range(num_scenarios):
            zetapbar=zetapbar_list[jsc]
            zetaplot=zetapbar[:,j].reshape((Npts1,Npts2),order='F')
            nanplot=isnan(zetaplot)
            notnanplot=~isnan(zetaplot)
            #diffplot=where(logical_and(nan0,nanplot),nan,0.0)
            #diffplot=where(logical_and(nan0,notnanplot),zetaplot,diffplot)
            diffplot=where(logical_and(nan0,notnanplot),zetaplot,0.0)
            diffplot=where(logical_and(notnan0,nanplot),zetaplot0,diffplot)
            diffplot=where(logical_and(notnan0,notnanplot),abs(zetaplot-zetaplot0),diffplot)

            #HERE
            ##  ensure zetadiffplot_weighted has no nans in it. If both are nans, no difference.
            ##  so using 0.
            zetadiffplot_weighted[jsc] += diffplot*(pbarlist[j]/sumpbarlist)
            ##  now add in the nans for plotting purposes, this is where did not get wet at all
            diffplot=where(logical_and(nan0,nanplot),nan,diffplot)
            #HERE

            #diffplot_pos = where(diffplot > 0.0,diffplot,nan)
            nanmax_zetadiffplot=nanmax(diffplot)
            #nanmin_zetadiffplot=nanmin(diffplot_pos)
            nanmean_zetadiffplot=nanmean(diffplot)

            masking=isnan(diffplot)
            diffplot=ma.masked_where(masking,diffplot)

            cmax = ceil(max(nanmax_zetadiffplot,5.)*100.)/100.

            outfile.write('\n')
            outfile.write('max value in zetadiffplot in %s was  %9.5f \n' \
                      %(legendnames[jsc],nanmax_zetadiffplot))
            #outfile.write('min value in positive zetadiffplot in %s was  %9.5f \n' \
            #         %(legendnames[jsc],nanmin_zetadiffplot))
            outfile.write('mean of zetadiffplot in %s was  %9.5f \n' \
                      %(legendnames[jsc],nanmean_zetadiffplot))

            if (cmax <= 7.0):
                clines = [1e-3] + [0.25,0.5,1.0,1.5] + list(linspace(2.0,7.0,6))
            else:
                clines = [1e-3] + [0.25,0.5,1.0,1.5] + list(linspace(2.0,6.0,5)) + [cmax]

            nlines = len(clines)
            n1 = int(floor((nlines-1)/2.))
            n2 = nlines - 1 - n1
            Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
            Red = hstack([linspace(0,0.8,n1), ones(n2)])
            Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
            colors2 = zip(Red,Green,Blue)

            #As jsc marches through the scenarios, the contour plots are plotted
            #using plots 2 columns wide, top to bottom, left to right
            plt_number=jsc+1
            plt.subplot(scenario_rows,2,plt_number) 
            CS100=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],diffplot,\
                       clines,colors=colors2)
            if (num_scenarios == 2):
                cbar_shrink=.7
            else:
                cbar_shrink=.9
            cbar = plt.colorbar(shrink=cbar_shrink)
            cbar.set_ticks(clines)
            plt.hold(True)

            # Contours of topo added:
            if 1:
                #plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
                plt.contour(X_CC,Y_CC,Z_CC, arange(0,8,4), colors='k')
                plt.ticklabel_format(format='plain',useOffset=False)
                #plt.xlim(-124.215,-124.18)
                plt.xlim(region_extent[0],region_extent[1])
                if (logical_or(plt_number == num_scenarios,plt_number == num_scenarios-1)):
                    plt.xticks(rotation=20,fontsize=10)
                    plt.yticks(fontsize=10)
                else:
                    plt.xticks([])
                    plt.yticks(fontsize=10)
                #plt.ylim(41.735,41.768)
                plt.ylim(region_extent[2],region_extent[3])

            if 0:
                #Draw green box around the xp,yp sub-contour plot
                box_x=[region_extent[0],region_extent[1],region_extent[1],\
                       region_extent[0],region_extent[0]]
                box_x=array(box_x)
                box_y=[region_extent[2],region_extent[2],region_extent[3],\
                      region_extent[3],region_extent[2]]
                box_y=array(box_y)
                plt.plot(box_x,box_y,'g-',linewidth='2')

            a = plt.gca()
            a.set_aspect(1./cos(aspectNO*pi/180.))
            titlestr='\nzeta-diff-contours, p= ' + str(pbarlist[j]) + ', ' + legendnames[jsc]
            plt.title(titlestr,fontsize=10)

        fname=Eventpath + '/zeta-diff-contours_prob_' + titlelist[j] + '.png'
        plt.savefig(fname)
        plt.clf()        #Clear figure 14 after all scenario plots are on the figure and saved

        todo=1
        if (logical_and(depth==True,todo==1)):
            for jsc in range(num_scenarios):
                zetapbar=zetapbar_list[jsc]
                zetaplot=zetapbar[:,j]
                zetaplot=zetapbar[:,j].reshape((Npts1,Npts2),order='F')
                B0_pos = B0_pos_list[jsc].reshape((Npts1,Npts2),order='F')

                #masking=logical_or(isnan(zetaplot),zetaplot <=0.0)
                masking=isnan(zetaplot)
                zetaplot=ma.masked_where(masking,zetaplot)
                Bplot=ma.masked_where(masking,B0_pos) #mask Bplot at same place as zetaplot
                                                      #Bplot can be masked differently each
                                                      #pass thru the loop, so call it Bplot
                                                      #Bplot will be masked 
                                                      #where no flood for this probability
                                                      #zeta=0.0 was not even EXCEEDED.
                cmax=nanmax(zetaplot+Bplot)
                #cmin=nanmin(zetaplot+Bplot)

                outfile.write('\n')
                outfile.write('max value in etaplot in %s was %9.5f \n' \
                      %(legendnames[jsc],cmax))
                if (cmax <= 16.0):
                    clines = [1e-3] + [.5,1.0] + list(linspace(2.0,16.0,8))
                else:
                    clines = [1e-3] + [.5,1.0] + list(linspace(2.0,14.0,7)) + [cmax]
                nlines = len(clines)
                n1 = int(floor((nlines-1)/2.))
                n2 = nlines - 1 - n1
                Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
                Red = hstack([linspace(0,0.8,n1), ones(n2)])
                Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
                colors2 = zip(Red,Green,Blue)

                plt_number=jsc+1
                plt.subplot(scenario_rows,2,plt_number) 
                CS100=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],zetaplot+Bplot,\
                            clines,colors=colors2)
                if (num_scenarios == 2):
                    cbar_shrink=.7
                else:
                    cbar_shrink=.9
                cbar = plt.colorbar(shrink=cbar_shrink)
                cbar.set_ticks(clines)
                plt.hold(True)

                # Contours of topo added:
                if 0:  #Would have to choose the B and shape before use here
                    clines2 = linspace(0,20,11)
                    plt.contour(xp[:-1,:-1],yp[:-1,:-1],B,clines2,colors='k')
                    plt.ticklabel_format(format='plain',useOffset=False)
                    plt.xticks(rotation=20)
                    plt.axis(region_extent)
                    a = plt.gca()
                    a.set_aspect(1./cos(aspectNO*pi/180.))
                if 1:
                    #plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
                    plt.contour(X_CC,Y_CC,Z_CC, arange(0,8,4), colors='k')
                    plt.ticklabel_format(format='plain',useOffset=False)
                    #plt.xlim(-124.215,-124.18)
                    plt.xlim(region_extent[0],region_extent[1])
                    if (logical_or(plt_number == num_scenarios,plt_number == num_scenarios-1)):
                        plt.xticks(rotation=20,fontsize=10)
                        plt.yticks(fontsize=10)
                    else:
                        plt.xticks([])
                        plt.yticks(fontsize=10)
                    #plt.ylim(41.735,41.768)
                    plt.ylim(region_extent[2],region_extent[3])

                if 0:
                    #Draw green box around the xp,yp contour plot
                    box_x=[region_extent[0],region_extent[1],region_extent[1],\
                           region_extent[0],region_extent[0]]
                    box_x=array(box_x)
                    box_y=[region_extent[2],region_extent[2],region_extent[3],\
                           region_extent[3],region_extent[2]]
                    box_y=array(box_y)
                    plt.plot(box_x,box_y,'g-',linewidth='2')

                a = plt.gca()
                a.set_aspect(1./cos(aspectNO*pi/180.))
                titlestr='\neta-contours for p= ' + str(pbarlist[j]) + ', '+legendnames[jsc]
                plt.title(titlestr,fontsize=10)
            fname=Eventpath + '/eta-contours_prob_' + titlelist[j] + '.png'
            plt.savefig(fname)
            plt.clf()   #clear figure 14

            ##Finding the difference of eta for all scenarios from eta of scenario jsc=0
            zetapbar0=zetapbar_list[0]
            zetaplot0=zetapbar0[:,j].reshape((Npts1,Npts2),order='F')
            B0_pos0 = B0_pos_list[0].reshape((Npts1,Npts2),order='F')
            etaplot0=zetaplot0+B0_pos0     #nothing masked yet, can add nans, or nan to a number
                                           #will be nan where zetaplot0 is nan. B0_pos0 never nan
            nan0=isnan(etaplot0)
            notnan0=~isnan(etaplot0)

            for jsc in range(num_scenarios):
                zetapbar=zetapbar_list[jsc]
                zetaplot=zetapbar[:,j]
                zetaplot=zetapbar[:,j].reshape((Npts1,Npts2),order='F')
                B0_pos = B0_pos_list[jsc].reshape((Npts1,Npts2),order='F')
                etaplot = zetaplot + B0_pos     #nothing masked yet
                nanplot = isnan(etaplot)
                notnanplot = ~isnan(etaplot)
                diffplot=where(logical_and(nan0,nanplot),nan,0.0)
                diffplot=where(logical_and(nan0,notnanplot),etaplot,diffplot)
                diffplot=where(logical_and(notnan0,nanplot),etaplot0,diffplot)
                diffplot=where(logical_and(notnan0,notnanplot),abs(etaplot-etaplot0),diffplot)
                #diffplot_pos=where(diffplot > 0.0,diffplot,nan)
                cmax=nanmax(diffplot)
                #cmin=nanmin(diffplot_pos)
                cave=nanmean(diffplot)
                masking=isnan(diffplot)
                diffplot=ma.masked_where(masking,diffplot)
                outfile.write('\n')
                outfile.write('max value in etadiffplot in %s was %9.5f \n' \
                      %(legendnames[jsc],cmax))
                #outfile.write('min value in positive etadiffplot in %s was %9.5f \n' \
                #     %(legendnames[jsc],cmin))
                outfile.write('mean of etadiffplot in %s was %9.5f \n' \
                      %(legendnames[jsc],cave))
                if (cmax <= 7.0):
                    clines = [1e-3] + [0.25,0.5,1.0,1.5] + list(linspace(2.0,7.0,6))
                else:
                    clines = [1e-3] + [0.25,0.5,1.0,1.5] + list(linspace(2.0,6.0,5)) + [cmax]

                nlines = len(clines)
                n1 = int(floor((nlines-1)/2.))
                n2 = nlines - 1 - n1
                Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
                Red = hstack([linspace(0,0.8,n1), ones(n2)])
                Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
                colors2 = zip(Red,Green,Blue)

                plt_number=jsc+1
                plt.subplot(scenario_rows,2,plt_number) 
                CS100=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],diffplot,\
                            clines,colors=colors2)
                if (num_scenarios == 2):
                    cbar_shrink=.7
                else:
                    cbar_shrink=.9
                cbar = plt.colorbar(shrink=cbar_shrink)
                cbar.set_ticks(clines)
                plt.hold(True)

                # Contours of topo added:
                if 1:
                    #plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
                    plt.contour(X_CC,Y_CC,Z_CC, arange(0,8,4), colors='k')
                    plt.ticklabel_format(format='plain',useOffset=False)
                    #plt.xlim(-124.215,-124.18)
                    plt.xlim(region_extent[0],region_extent[1])
                    if (logical_or(plt_number == num_scenarios,plt_number == num_scenarios-1)):
                        plt.xticks(rotation=20,fontsize=10)
                        plt.yticks(fontsize=10)
                    else:
                        plt.xticks([])
                        plt.yticks(fontsize=10)
                    #plt.ylim(41.735,41.768)
                    plt.ylim(region_extent[2],region_extent[3])

                if 0:
                    #Draw green box around the xp,yp contour plot
                    box_x=[region_extent[0],region_extent[1],region_extent[1],\
                           region_extent[0],region_extent[0]]
                    box_x=array(box_x)
                    box_y=[region_extent[2],region_extent[2],region_extent[3],\
                           region_extent[3],region_extent[2]]
                    box_y=array(box_y)
                    plt.plot(box_x,box_y,'g-',linewidth='2')

                a = plt.gca()
                a.set_aspect(1./cos(aspectNO*pi/180.))
                titlestr='\neta-diff-contours for p= ' + str(pbarlist[j]) + ', '+legendnames[jsc]
                plt.title(titlestr,fontsize=10)
            fname=Eventpath + '/eta-diff-contours_prob_' + titlelist[j] + '.png'
            plt.savefig(fname)
            plt.clf()   #clear figure 14
        #if logical or is now done

    ##HERE, at this point zetadiffplot_weighted is done for all the scenarios. We can use figure 14
    ##      to plot, or we could ask for the expected value (weighted value) at particular points to
    ##      be printed out for each scenario, or we could do both.
    for jsc in range(num_scenarios):
        maxvalue=zetadiffplot_weighted[jsc].max()
        cmax = ceil(max(maxvalue,4.))*100./100.
        meanvalue=zetadiffplot_weighted[jsc].mean()
        outfile.write('\n')
        outfile.write('max value in zetadiffplot_weighted in %s was  %9.5f \n' \
                      %(legendnames[jsc],maxvalue))
        outfile.write('mean of zetadiffplot_weighted in %s was  %9.5f \n' \
                      %(legendnames[jsc],meanvalue))
        if (cmax <= 4.0):
            clines = [1e-3] + [0.1,0.25,0.5,1.0] + list(linspace(1.5,4.0,6))
        else:
            clines = [1e-3] + [0.1,0.25,0.5,1.0] + list(linspace(1.5,3.5,5)) + [cmax]

        nlines = len(clines)
        n1 = int(floor((nlines-1)/2.))
        n2 = nlines - 1 - n1
        Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
        Red = hstack([linspace(0,0.8,n1), ones(n2)])
        Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
        colors2 = zip(Red,Green,Blue)

        #As jsc marches through the scenarios, the contour plots are plotted
        #using plots 2 columns wide, top to bottom, left to right
        plt_number=jsc+1
        plt.subplot(scenario_rows,2,plt_number) 
        CS100=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],zetadiffplot_weighted[jsc],\
                   clines,colors=colors2)
        if (num_scenarios == 2):
            cbar_shrink=.7
        else:
            cbar_shrink=.9
        cbar = plt.colorbar(shrink=cbar_shrink)
        cbar.set_ticks(clines)
        plt.hold(True)

        # Contours of topo added:
        if 1:
            #plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
            plt.contour(X_CC,Y_CC,Z_CC, arange(0,8,4), colors='k')
            plt.ticklabel_format(format='plain',useOffset=False)
            #plt.xlim(-124.215,-124.18)
            plt.xlim(region_extent[0],region_extent[1])
            if (logical_or(plt_number == num_scenarios,plt_number == num_scenarios-1)):
                plt.xticks(rotation=20,fontsize=10)
                plt.yticks(fontsize=10)
            else:
                plt.xticks([])
                plt.yticks(fontsize=10)
            #plt.ylim(41.735,41.768)
            plt.ylim(region_extent[2],region_extent[3])

        if 0:
            #Draw green box around the xp,yp sub-contour plot
            box_x=[region_extent[0],region_extent[1],region_extent[1],\
                   region_extent[0],region_extent[0]]
            box_x=array(box_x)
            box_y=[region_extent[2],region_extent[2],region_extent[3],\
                  region_extent[3],region_extent[2]]
            box_y=array(box_y)
            plt.plot(box_x,box_y,'g-',linewidth='2')

        a = plt.gca()
        a.set_aspect(1./cos(aspectNO*pi/180.))
        titlestr='\nzeta-diff-contours-weighted ' + ', ' + legendnames[jsc]
        plt.title(titlestr,fontsize=10)
    fname=Eventpath + '/zeta-diff-contours-weighted.png'
    plt.savefig(fname)
    plt.clf()        #Clear figure 14 after all scenario plots are on the figure and saved



    #NOW print the value of zetadiffplot_weighted for the rows in whichrows for
    #    each scenario
    outfile.write('\n')
    for irow in whichrows:
        outfile.write(' FOR IROW=: %s ' %irow)
        for jsc in range(num_scenarios):
            zetadiffplot_weighted_linear=zetadiffplot_weighted[jsc].reshape(Npts1*Npts2,order='F')
            thevalue=zetadiffplot_weighted_linear[irow]
            outfile.write('scenario: %s zetadiffplot_weighted: %9.2f \n' %(legendnames[jsc],thevalue))
        outfile.write('\n')
    #end of printing of zetadiffplot_weighted, for each row used for the Hazard Curves

    plt.close(14)
    outfile.close()
    #raw_input(' hit return to end session ')
    return zetapbar

    
if __name__ == '__main__':
###### TYPICAL DATA  #########################################################

    #comparison_number=1
    #scenarionames=['scenario_all_fine_uniform','scenario_uniform10']
    #scenarioresolution=['fine','fine']

    #comparison_number=2
    #scenarionames=['scenario_all_fine_uniform','scenario_uniform10',\
    #              'scenario_all_fine','scenario_highest10']
    #scenarioresolution=['fine','fine','fine','fine']

    #comparison_number=3
    #scenarionames=['scenario_all_fine','scenario_highest10']
    #scenarioresolution=['fine','fine']

    comparison_number=4
    scenarionames=['scenario_all_fine','scenario_all_coarse']
    scenarioresolution=['fine','coarse']
################################################################################
    projectdir = os.environ['FEMA']

    print 'INSIDE name is main DRIVER, data echo: '
    print ' '
    print 'The project directory was: ',projectdir
    print ' '
    print 'The scenarionames were: '
    print scenarionames
    print ' '
    print 'The comparison number for this comparison was: '
    print comparison_number
    print ' '

    #Read in the file below only once, since same one used for all transects and 2Dgrids
    fixed_big_grid_file = projectdir + '/DataFiles/' + 'CCtopo_1_3sec.tt3'
    print 'The fixed grid file for the biggest fixed grid region used for plotting '
    print 'little map in upper right graph corner: '
    print fixed_big_grid_file
    print ' '
    from clawpack.geoclaw.topotools import Topography
    CCtopo = Topography()
    CCtopo.read(fixed_big_grid_file,3)

    scenariodir=[]               #going to be a list
    fg_types_todo=[]             #going to be a list of arrays
    j=0
    for scenarioname in scenarionames:
        scenariodir.append(projectdir + '/' + scenarioname)
        fg_types_todo_file = projectdir + '/' + scenarioname + '/fg_types_todo.txt'
        types_array=loadtxt(fg_types_todo_file,dtype=integer)
        fg_types_todo.append(types_array)
        j+=1

    print 'the scenariodir file was: '
    print scenariodir
    print ' '
    print 'the fg_types_todo list was: '
    print fg_types_todo
    print ' '

    fg_2Dgrids_todo_file=[]
    num_scenarios = len(scenarionames)
    for j in range(num_scenarios):
        if (set(list(fg_types_todo[0])) == set(list(fg_types_todo[j]))): #scenarios same types
            for itype in fg_types_todo[0]:
                if (itype == 2):
                    fg_2Dgrids_todo_file.append(scenariodir[j] + '/fg_2Dgrids_todo.txt')
        else:
            print 'can not compare these scenarios. Their todo files have different types'
            print 'scenario number j different type from scenario number 0, j= ',j
            print 'fg_types_todo[0] was: ',fg_types_todo[0]
            print 'fg_types_todo[j] was: ',fg_types_todo[j]
            sys.exit(1)

    print 'the fg_2Dgrids_todo_file list was: '
    print fg_2Dgrids_todo_file
    print ' '

    ######  Now loop over all the work that has been done for these scenarios, and ###
    ######  check which 2Dgrids to compare                                         ###
    found=0
    for itype in fg_types_todo[0]:    #set(list(fg_types_todo[j]))=set(list(fg_types_todo[0]))##
        if (itype == 2):              #Do all the 2Dgrid type of fixed grids
            found=1
            fg_2Dgrids_todo=[]        #A list of arrays
            for j in range(num_scenarios):
                fg_2Dgrids_todo.append(loadtxt(fg_2Dgrids_todo_file[j],dtype=integer))
                if (size(fg_2Dgrids_todo[j]) == 1):
                    fg_2Dgrids_todo[j] = array([fg_2Dgrids_todo[j]])
                if (set(list(fg_2Dgrids_todo[0])) != set(list(fg_2Dgrids_todo[j]))):
                    print 'can not compare 2Dgrids of these scenarios. '
                    print 'because they do not have the same 2Dgrid numbers'
                    print 'fg_2Dgrids_todo[0] was: ',fg_2Dgrids_todo[0]
                    print 'j and fg_2Dgrids_todo[j] were: ',j,fg_2Dgrids_todo[j]
                    sys.exit(1)

            #All 2Dgrids have the same numbers, so continue
            print 'the fg_2Dgrids_todo[0] was: '
            print fg_2Dgrids_todo[0]
            print ' '

            for iwhich in fg_2Dgrids_todo[0]:       #all fg_2Dgrids_todo are equal as a set
                print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                print '++++++++++++ PROCESSING 2Dgrid number (iwhich) +++++++ ',iwhich
                print ' '
                hazard_curves_file=[]; OutZetafile=[]; #2Dgrid iwhich, multiple scenarios
                fixed_grid_file=[]
                for j in range(num_scenarios):
                    outdir=scenariodir[j]+'/'+scenarionames[j]+'_2Dgrid_fg'+str(iwhich)
                    hazard_curves_file.append(outdir + '/hazardcurves_probs.npy')
                    OutZetafile.append(outdir + '/' + 'ZetaFloods.npy')
                    print 'The hazard curves .npy file for jth scenario when j= ',j,' was: '
                    print hazard_curves_file[j]
                    print ' '
                    print 'The Flood zetapbar .npy file for jth scenario when j= ',j,' was: '
                    print OutZetafile[j]
                    print ' '
                    type_name=scenarioresolution[j]
                    fg_file = projectdir+'/' + type_name +'_runs/'+\
                                       type_name+'_xyB0_fg'+str(iwhich) +'.npy'
                    fixed_grid_file.append(fg_file)
                    print 'The fixed_grid_file .npy file for the jth scenario when j= ',j,' was: '
                    print fixed_grid_file[j]
                    print ' '

                ###Note, need to read Npts1 and Npts2 from some file region_dimensions.txt
                ###      in outdir that make_scenerio.py put there.
                ###      The file has length two with the values  Npts1 and  Npts2
                ###
                ###All the scenarios will have the same Npts1 and Npts2, so just retrieve
                ###this information from the first one.

                outdir=scenariodir[0]+'/'+scenarionames[0]+'_2Dgrid_fg'+str(iwhich)
                dimension_file=outdir + '/' + 'region_dimensions.txt'
                dimensions=loadtxt(dimension_file,dtype=integer)
                Npts1=dimensions[0]; Npts2=dimensions[1];
                print 'The 2Dgrid was Npts1 x Npts2 points in longxlat dimensions'
                print Npts1,Npts2
                print ' '

                #The plotter will automatically choose evenly distributed points for
                #the hazard curves.  If these are good enough, choose whichrows=None;
                #otherwise, make whichrows a list (say of gauge indices) and the plotter
                #will append the evenly distributed points to the whichrows list.
                #BUT NOTE THAT the gauge indices can be 2Dgrid dependent, and if this
                #is a priority, we would read them in from right scenario directory,
                #and make_scenario.py would have to be modified to put them there.

                whichrows=None
                depth=True; speed=False; mflux=False;
                zeta = linspace(0,12,1201)
                whichzetas=[0,100,200,300,400,600,800,1000]  #zeta=0.0,1.0,2.0,3.0,4.0,6.0,8.0,10.0
                whichzetas=array(whichzetas)
                pbarlist=[.1,.2,.3,.4,.5,.6,.7,.8,.9,.95]   #10 zeta-contours 

                comparisondir=projectdir + '/comparisons/comparison_' + str(comparison_number) +\
                              '/2Dgrid_fg' + str(iwhich)
                print 'The comparison outputs for this 2Dgrid will be written at: '
                print comparisondir
                print ' '
                print 'CALLING hazard_curves_2Dgrid_compare for fg 2Dgrid: ',iwhich
                print ' '

                if 1:
                    hazard_curves_2Dgrid_compare(hazard_curves_file,fixed_grid_file,\
                          OutZetafile,\
                          zeta,whichzetas,comparisondir,pbarlist,Npts1,Npts2,\
                          CCtopo.X,CCtopo.Y,CCtopo.Z,scenarionames,scenarioresolution,\
                          whichrows=whichrows,speed=speed,depth=depth,mflux=mflux)
    if (found == 0):
        print 'No 2Dgrids were found to compare'
