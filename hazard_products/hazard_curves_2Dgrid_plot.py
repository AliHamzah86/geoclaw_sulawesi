"""

Usage: 
 >>>
 from hazard_curves_2Dgrid_plot import hazard_curves_2Dgrid_plot

    where for example

    INPUT: 
      hazard_curves_file = the .npy file used for the plots:

           P1 P2 P3 PNzeta 

      fixed_grid_file = file with x y B in colms, a .npy file
        (for the FEMA project, B is B0, initial bathy)

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


      The following region_extent will be calculated in the routine.

      region_extent: [longwest,longeast,latsouth,latnorth] for the 2D fixed grid. Use
                     for plotting.


      Example:
         import os
         CCdir = os.environ['CC']
         probsnpy='hazardcurves_probs.npy'
         Testsdir=CCdir + '/Realizations_Test/'
         fixed_grid_file = CCdir + '/Datafiles/fixedgrid_xyB_test.npy'

         hazard_curves_file  = Testsdir +'EVENTeall' + '/'+ probsnpy 

         See the __main__ below for setting zeta, whichrows, Npts1, Npts2
         and region_extent

         whichrows is now set within the routine below based on the number
                   of grid points on the fixed grid - done automatically.
                   Particular gauges could be added in the calling program
                   by setting whichrows as a list with these gauge indices,
                   then this program appends to the inputted list, the
                   automatic numbers.


         speed=True if we are doing plots of speeds; otherwise speed=False
         depth=True if we are doing plots of depth; otherwise depth=False
         mflux=True if we are doing plots of mflux; otherwise mflux=False
   

    OUTPUT:  Does a plot for each row of the hazard_curves_file that
             is indicated in whichrows.  The first row is row 0.
             User may be prompted to step through the plots.

             Also does zeta-contour contourf plots
             where p is specified.  Say the p=.01 "event" (100 year flood)
             Also could writes the  Zeta values for these two probs
             on the files below for further examination.

             For example:

             outZetafile=Testsdir + 'EVENTeall' + '/' + 'Zeta100_500.npy'

             Also does p-contour  plots of the probabilities
             associated with a particular zetai level as specified in
             the list whichzetas. User may be prompted to step through the
             plots. The list whichzetas are the columns in the 
             hazard_curves_file that we want: 

             For example:

             whichzetas=[0,5]

             would correspond to zeta0 and zeta5. 


"""
from numpy import loadtxt
from numpy import array,shape,load,save,nan,isnan,linspace,sum,where
from numpy import zeros,ones,ma,amin,amax,logical_or,logical_and
from numpy import int,integer,size,floor,ceil
from numpy import cos,pi,hstack,arange
from matplotlib.mlab import find
from matplotlib import image,colors
import matplotlib.pyplot as plt
import os,sys

def hazard_curves_2Dgrid_plot(hazard_curves_file,fixed_grid_file,\
       zeta,whichzetas,Eventpath,pbarlist,Npts1,Npts2,\
       X_CC,Y_CC,Z_CC,\
       whichrows=None,speed=False,depth=True,mflux=False):

    #Specific to Crescent City and the FEMA project
    #FEMA=os.environ['FEMA']                           #If we need this
    #CCmap=image.imread(FEMA+'/maps/CCmap8.png')       #don't have this yet

    aspectNO=41.75                                     #aspect latitude for plotting
                                                       #specific to Crescent City

    commentsfile=Eventpath + '/hc_plot_run.output'
    outfile=open(commentsfile,'w')

    outfile.write(' THE 2D fixed grid dimensions Npts1 and Npts2 are = : ')
    outfile.write('%s %s ' %(Npts1,Npts2))
    outfile.write(' \n')
    outfile.write(' \n')
    #
    fprobs=load(hazard_curves_file)
    d=load(fixed_grid_file)
    x=d[:,0]
    y=d[:,1]
    B=d[:,2]           #This is B0 for the FEMA project

    #Compute the region_extent from the x and y information
    region_extent=[]; region_extent.append(x[0]); region_extent.append(x[-1]);
    region_extent.append(y[0]); region_extent.append(y[-1]);
    outfile.write(' THE region_extent for this fgno fixed grid is : \n')
    outfile.write('%s ' %region_extent)
    outfile.write(' \n')
    outfile.write(' \n')

    if (depth==True):
        outfile.write(' This is a DEPTH zeta run.  \n')
        outfile.write(' Plots in directory: %s ' %Eventpath)
        outfile.write(' \n')
        outfile.write(' \n')
    elif (speed==True):
        outfile.write(' This is a SPEED zeta run.  \n')
        outfile.write(' Plots in directory: %s ' %Eventpath)
        outfile.write(' \n')
        outfile.write(' \n')
    elif (mflux==True):
        outfile.write(' This is a MFLUX zeta run.  \n')
        outfile.write(' Plots in directory: %s ' %Eventpath)
        outfile.write(' \n')
        outfile.write(' \n')
    else:
        print ('program stopping, depth,speed,mflux all False')
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

    #getting B ready for plotting
    #getting the x and y ready for pcolor or contourf plots as xp and yp
    #Note: We only need to do this is the fixed grid is of type fgno=3;
    #      otherwise, we use B and y in the plotting, no need for xp at all. 
    #
    B=B.reshape((Npts1,Npts2),order='F')
    dx=.000001/2.; dy=dx;
    xp=zeros((Npts1+1,Npts2+1))
    yp=zeros((Npts1+1,Npts2+1))
    xp[0:Npts1,0:Npts2]=x.reshape((Npts1,Npts2),order='F')
    xp[0:Npts1,Npts2]=xp[0:Npts1,Npts2-1]        #add another colm for same x-value
    xp[Npts1,:]=xp[Npts1-1,:]+dx                 #add another row with bigger x-value
    yp[0:Npts1,0:Npts2]=y.reshape((Npts1,Npts2),order='F')
    yp[0:Npts1,Npts2]=yp[0:Npts1,Npts2-1]+dy     #add another col with bigger y-value
    yp[Npts1,:]=yp[Npts1-1,:]                    #add another row for same y-value           

    plt.figure(11,(12,8))
    fprobs121=zeros(121)        
    for irow in whichrows:       #Grid Point Hazard Curve Plots
        plt.clf()
        fprobs121=fprobs[irow,:]
        masking=where(isnan(fprobs121),True,False)
        fprobs121=ma.masked_where(masking,fprobs121)


        #fprobs121=where(fprobs121 <= 1e-10,1e-10,fprobs121) #So semilog can work

        if (mflux == True):
            fprobs121=where(fprobs121 <= 1e-10,1e-10,fprobs121) #So semilog can work
            plt.loglog(zeta,fprobs121,'ro')
            plt.hold(True) 
            plt.loglog(zeta,fprobs121,'r-')
        else:                                  #Not doing semilog for depth or speed plots 
                                               #We don't have recurrence time of earthquake
                                               #We just have conditional probs of realizations
            #plt.plot(zeta,fprobs121,'ro',markersize=5)
            #plt.hold(True) 

            plt.plot(zeta,fprobs121,'r-')
            #plt.semilogy(zeta,fprobs121,'ro')
            #plt.hold(True) 
            #plt.semilogy(zeta,fprobs121,'r-')

        #ylim(0,.040)
        #ylim(0,.010)
        plt.xlim(0,zeta[-1])
        plt.ylim(-0.1,1.1) 

        #Not using semilogy so comment this out and choose ylim automatically
        #plt.ylim(1e-10,100)

        #Use yticks for semilogy plots
        #yticks([1e-5,1e-4,1e-3,.002,.01,.1,1],\
        #       ['1e-5','1e-4','1e-3','0.002','0.01','0.1','1'])
        #yticks([1e-10,1e-7,1e-6,1e-5,1e-4,.0004,1e-3,.002,.01,.1,1],\
        #  ['1e-10','1e-7','1e-6','1e-5','1e-4','.0004','1e-3','0.002','0.01','0.1','1'])

        #Use yticks for tideprobs, or for the conditional probs of realizations
        plt.yticks([0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0],\
           ['0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'])
        plt.xticks([0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0],\
           ['0','2','4','6','8','10','12','14','16','18','20','22','24'])

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
        plt.title('At Longitude, Latitude, Bathy:  (%9.4f,%9.4f,%9.4f) ' \
                                                      %(x[irow],y[irow],d[irow,2]))

        ######### putting little map on the plot, showing the location
        plt.axes([.7,.65,.25,.25]) # .7,.65 is lower left proportion of figure, .25,.25 width and ht
                                   # proportion of figure
        plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
        a = plt.gca()
        a.set_aspect(1./cos(aspectNO*pi/180.))

        xi=x[irow]; yi=y[irow];
        plt.plot([xi],[yi],'ro')
        plt.xticks([])
        plt.yticks([])
        #########


        #########

        #raw_input(' hit return to get the next plot ')
        plt.savefig(fname)

    plt.figure() 
    for izeta in whichzetas:      #p-contour plots for heights in whichzetas
        plt.clf()
        zetarest=fprobs[:,izeta].reshape((Npts1,Npts2),order='F')
        masking=logical_or(isnan(zetarest),zetarest <=0.0)
        zetarest=ma.masked_where(masking,zetarest)
        #
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
        cbar = plt.colorbar()
        cbar.set_ticks(V)
        plt.hold(True)

        # Contours of topo added:
        if 0:
            clines2 = linspace(0,20,11)
            plt.contour(xp[:-1,:-1],yp[:-1,:-1],B,clines2,colors='k')
            plt.ticklabel_format(format='plain',useOffset=False)
            plt.xticks(rotation=20)
            plt.axis(region_extent)
            a = plt.gca()
            a.set_aspect(1./cos(aspectNO*pi/180.))
        if 1:
            plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
            plt.xticks(rotation=20)
            plt.xlim(-124.215,-124.18)
            plt.ylim(41.735,41.768)
            plt.ticklabel_format(format='plain',useOffset=False)

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

        nn=str(int(zeta[izeta]*100)).zfill(5)
        if (depth==True):
            plt.title('p-contours for zeta=%5.2f m' %(zeta[izeta]), fontsize=20)
            fname=Eventpath + '/p-contours_zeta_%s.png' %nn
        elif (speed==True):
            plt.title('p-contours for speed=%5.2f m/sec' %(zeta[izeta]), fontsize=20)
            fname=Eventpath + '/p-contours_speed_%s.png' %nn
        else:
            plt.title('p-contours for momentum flux=%5.2f ' %(zeta[izeta]), fontsize=20)
            fname=Eventpath + '/p-contours_mflux_%s.png' %nn
        plt.savefig(fname)

        if 0:   #Don't do over a map at the moment
            #raw_input(' hit return for p-contourf over map for zeta= %5.2f ' %(zeta[izeta]))
            plt.clf()
            #extent = (360-124.23, 360-124.16, 41.73,41.785)
            #extent = (235.77, 235.85, 41.7348,41.785)

            #WOULD NEED A BETTER MAP OF OUR FEMA REGION to use statement below
            plt.imshow(CCmap,extent=region_extent) 

            #V=[0.0,.0004,.001,.0015,.002,.0025,.005,.0075,.01,.015,.02,1.]
            #V=[0.0,.0004,.001,.002,.004,.005,.01,.02,.04,.06,1.] used for CCdir

            #Use for tideprobs for Discovery and Lynch or conditional prob situations like FEMA
            V=[0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]
            V2=[0.0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.000001]
            #V=[0.0,.0004,.001,.002,.005,.01,.02,.04,.1,.2,1.]
            #V=[0.0,.00001,.0001,.0002,.0004,.0008,.0010,.0012,.0014,.0016,.0018,.0020,.0040]
            norm=colors.BoundaryNorm(V,256)
            CS=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],zetarest,V2,norm=norm,alpha=0.8)
            cbar = plt.colorbar()
            cbar.set_ticks(V)
            plt.contour(xp[:-1,:-1],yp[:-1,:-1],B,[0],colors='w')
            plt.contour(xp[:-1,:-1],yp[:-1,:-1],zetarest,[0.002,0.01],colors='k')
            plt.ticklabel_format(format='plain',useOffset=False)
            plt.xticks(rotation=20)
            plt.axis(region_extent)
            a = plt.gca()
            a.set_aspect(1./cos(aspectNO*pi/180.))
            nn=str(int(zeta[izeta]*100)).zfill(5)
            if (depth==True):
                plt.title('p-contours for zeta=%5.2f meters' %(zeta[izeta]), fontsize=20)
                fname=Eventpath + '/p-contours_map_zeta_%s.png' %nn
            elif (speed==True):
                plt.title('p-contours for speed=%5.2f m/sec' %(zeta[izeta]), fontsize=20)
                fname=Eventpath + '/p-contours_map_speed_%s.png' %nn
            else:
                plt.title('p-contours for momentum flux =%5.2f ' %(zeta[izeta]), fontsize=20)
                fname=Eventpath + '/p-contours_map_mflux_%s.png' %nn
            plt.savefig(fname)

    #Determine the zetapt01 and zetapt002 for the 100 yr and 500 yr events 
    #and store in the two columns of zetapbar
    #Find zetapt01 and zetapt002, zeta300, zeta975, zeta2500, zeta5000

    #Note, made the program more general, the particular "floods" are in pbarlist below

    sh=shape(fprobs)
    Npts=sh[0]
    Nzeta=sh[1]

    numpbar=len(pbarlist)
    zetapbar=nan*ones((Npts,numpbar))

    outfile.write('We compute fprobs-pbar for all Nzeta colms in fprobs \n')
    outfile.write('for each grid point. \n')

    okval=[]
    titlelist=[]

    for j in range(len(pbarlist)):
        titlelist.append('p=' + str(pbarlist[j]))
    
    j=0
    for pbar in pbarlist:
        outfile.write('  \n')
        outfile.write('pbar is: %10.6f' %pbar)
        outfile.write('  \n')
        fprobs_mpbar=fprobs-pbar

        num_nan=sum(where(isnan(fprobs_mpbar),1,0),1)
        num_pos=sum(where(fprobs_mpbar > 0.0,1,0),1)
        num_eq=sum(where(fprobs_mpbar == 0.0,1,0),1)
        num_neg=Nzeta - (num_pos + num_eq +num_nan)

        Ngridpts=Npts1*Npts2
        outfile.write('The number of grid points was: %5d \n' %Ngridpts)
        outfile.write('  \n')

        ind_all_nan=find(num_nan == Nzeta)
        outfile.write('no. of grid pts having fprobs-pbar=nan \
               in all colms is: %5d \n' %len(ind_all_nan))

        ind0=find(num_eq > 0.0)
        outfile.write('no. of grid pts having fprobs-pbar=0 \
               in some colm is: %5d \n' %len(ind0))
        zetapbar[ind0,j]=zeta[num_pos[ind0]+num_eq[ind0] -1]

        ind1=find(num_pos == Nzeta)
        outfile.write('no. of grid pts having fprobs-pbar>0 \
               in all colms is: %5d \n' %len(ind1))
        if (len(ind1) > 0):
            outfile.write('   BEWARE, zeta[-1] not big enough \n')
        zetapbar[ind1,j]=zeta[-1]

        ind2=find(num_neg == Nzeta)
        outfile.write('no. of grid pts having fprobs-pbar<0 \
               in all colms is: %5d \n' %len(ind2))
        zetapbar[ind2,j]=nan

        if (len(ind2)!= Npts): #checking this pbar is achieved
           okval.append(True)
        else:
           okval.append(False)

        ind3=find((num_eq == 0) & (num_pos > 0) & (num_neg > 0))
        outfile.write('no. of grid pts having no =, >=1 +  \
               and >=1 -: %5d \n' %len(ind3))
        left=num_pos[ind3]-1
        right=left+1
        #use fprobs[ind3,left] and fprobs[ind3,right]
        m=(zeta[right]-zeta[left])/(fprobs[ind3,right]-fprobs[ind3,left])
        zetapbar[ind3,j]=zeta[right]+m*(pbar-fprobs[ind3,right])
        j += 1

    outfile.write(' \n')
    outfile.write(' \n')
    outfile.write(' okval list was: %s ' %okval)
    outfile.write(' \n')
    outfile.write(' \n')

    #Note: Updated the program so pbar can be any value from
    #      the list pbarlist.
    #      okval is a list of True or False for whether there is
    #      a hazard flood for this value or not.  
    #      So need to check the okval list to know.
    #
    #
    # Set up and all the floods in the pbarlist if desired
    #
    B0_pos = where(B < 0.0,0,B)   #used to make eta. If a water point, B already was added in the
                                  #hmax used to make zetapbar, so zetaplot already has B effect, but
                                  #for land points we need to add B to zetaplot to get the eta.

    for j in range(len(pbarlist)):
        plt.clf()
        zetaplot=zetapbar[:,j].reshape((Npts1,Npts2),order='F')
        #masking=logical_or(isnan(zetaplot),zetaplot <=0.0)
        masking=isnan(zetaplot)   #will be consistent with the compare program
                                  #masked where "no flood", see compare program
        zetaplot=ma.masked_where(masking,zetaplot)
        Bplot=ma.masked_where(masking,B0_pos)      #masking will change each loop pass, so call
                                                   #the file to plot Bplot.
                                                   #Both zetaplot and Bplot are masked where there
                                                   #is no flood for this probability pbar

        if (okval[j] == True):
            cmax = ceil(max(zetaplot.max(), 5.)*100.)/100.
            outfile.write('\n')
            outfile.write(' FOR PBAR == : %s ' %pbarlist[j])
            outfile.write('max value in zetaplot and etaplot %7.3f %7.3f \n' \
                      %(amax(zetaplot),amax(zetaplot+Bplot)))
            outfile.write('min value in zetaplot and etaplot %7.3f %7.3f \n' \
                      %(amin(zetaplot),amin(zetaplot+Bplot)))
            outfile.write('\n')
            outfile.write('\n')

        else:
            cmax=5.
            outfile.write('\n')
            outfile.write(' FOR PBAR == : %s ' %pbarlist[j])
            outfile.write('No contour exists for this PBAR  \n')
            outfile.write('\n')
            outfile.write('\n')

        if 1:
            #raw_input(' hit return for this PBAR contourf plot ')
            #cmax = ceil(max(zetaplot.max(), 5.)*100.)/100.

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

                clines=[0.01, 0.05]+list(linspace(0.1,2,20))+[4,6,8,10,20,50,100,200,400,600,1000,dmax]
            nlines = len(clines)
            n1 = int(floor((nlines-1)/2.))
            n2 = nlines - 1 - n1
            Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
            Red = hstack([linspace(0,0.8,n1), ones(n2)])
            Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
            colors2 = zip(Red,Green,Blue)
            CS100=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],zetaplot,\
                       clines,colors=colors2)
            cbar=plt.colorbar()
            cbar.set_ticks(clines)
            plt.hold(True)

            # Contours of topo added:
            if 0:
                clines2 = linspace(0,20,11)
                plt.contour(xp[:-1,:-1],yp[:-1,:-1],B,clines2,colors='k')
                plt.ticklabel_format(format='plain',useOffset=False)
                plt.xticks(rotation=20)
                plt.axis(region_extent)
                a = plt.gca()
                a.set_aspect(1./cos(aspectNO*pi/180.))
            if 1:
                plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
                plt.ticklabel_format(format='plain',useOffset=False)
                plt.xticks(rotation=20)
                plt.xlim(-124.215,-124.18)
                plt.ylim(41.735,41.768)

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

            if (depth==True):
                titlestr='zeta-contours for p= ' + str(pbarlist[j])
                plt.title(titlestr)
                fname=Eventpath + '/zeta-contours_prob_' + titlelist[j] + '.png'
            elif (speed==True):
                titlestr='speed-contours for p= ' + str(pbarlist[j])
                plt.title(titlestr)
                fname=Eventpath + '/speed-contours_prob_' + titlelist[j] + '.png'
                #fname=Eventpath + '/speed-contours_prob_010.png'
            else:
                titlestr='mflux-contours for p= ' + str(pbarlist[j])
                plt.title(titlestr)
                fname=Eventpath + '/mflux-contours_prob_' + titlelist[j] + '.png'
                #fname=Eventpath + '/mflux-contours_prob_010.png'
            plt.savefig(fname)

        if 0:   #Not doing the plot over a map now
            #raw_input(' hit return for contourf plot over map for this PBAR ')
            plt.clf()
            #extent = (360-124.23, 360-124.16, 41.73,41.785)
            #extent = (235.77, 235.85, 41.7348,41.785)
            plt.imshow(CCmap,extent=region_extent)
            plt.contourf(xp[:-1,:-1],yp[:-1,:-1],zetaplot,clines,\
                     colors=colors2,alpha=0.8)
            cbar = plt.colorbar()
            cbar.set_ticks(clines)
            #clines3 = [1e-3,1,2]
            #plt.contour(xp[:-1,:-1],yp[:-1,:-1],zetaplot,clines3,colors='w')
            plt.contour(xp[:-1,:-1],yp[:-1,:-1],B,[0],colors='w')
            plt.ticklabel_format(format='plain',useOffset=False)
            plt.xticks(rotation=20)
            plt.axis(region_extent)
            a = plt.gca()
            a.set_aspect(1./cos(aspectNO*pi/180.))

            if (depth==True):
                titlestr='zeta-contours for p= ' + str(pbarlist[j])
                plt.title(titlestr)
                fname=Eventpath + '/zeta-contours_prob_' + titlelist[j] + '_map.png'
                #fname=Eventpath + '/zeta-contours_prob_010_map.png'
            elif (speed==True):
                titlestr='speed-contours for p= ' + str(pbarlist[j])
                plt.title(titlestr)
                fname=Eventpath + '/speed-contours_prob_' + titlelist[j] + '_map.png'
            else:
                titlestr='mflux-contours for p= ' + str(pbarlist[j])
                plt.title(titlestr)
                fname=Eventpath + '/mflux-contours_prob_' + titlelist[j] + '_map.png'
            plt.savefig(fname)

        todo=1
        if (logical_and(depth==True,todo==1)):
            plt.clf()
            cmax=(zetaplot+Bplot).max()

            if (cmax > 16.0):
                clines = [1e-3] + [.5,1.0] + list(linspace(2.0,14.0,7)) + [cmax]
            else:
                clines = [1e-3] + [.5,1.0] + list(linspace(2.0,16.0,8))
            nlines = len(clines)
            n1 = int(floor((nlines-1)/2.))
            n2 = nlines - 1 - n1
            Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
            Red = hstack([linspace(0,0.8,n1), ones(n2)])
            Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
            colors2 = zip(Red,Green,Blue)
            CS100=plt.contourf(xp[:-1,:-1],yp[:-1,:-1],zetaplot+Bplot,\
                       clines,colors=colors2)
            cbar=plt.colorbar()
            cbar.set_ticks(clines)
            plt.hold(True)

            # Contours of topo added:
            if 0:
                clines2 = linspace(0,20,11)
                plt.contour(xp[:-1,:-1],yp[:-1,:-1],B,clines2,colors='k')
                plt.ticklabel_format(format='plain',useOffset=False)
                plt.xticks(rotation=20)
                plt.axis(region_extent)
                a = plt.gca()
                a.set_aspect(1./cos(aspectNO*pi/180.))
            if 1:
                plt.contour(X_CC,Y_CC,Z_CC, arange(0,21,2), colors='k')
                plt.ticklabel_format(format='plain',useOffset=False)
                plt.xticks(rotation=20)
                plt.xlim(-124.215,-124.18)
                plt.ylim(41.735,41.768)

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

            titlestr='eta-contours (h+B0) for p= ' + str(pbarlist[j])
            plt.title(titlestr)
            fname=Eventpath + '/eta-contours_prob_' + titlelist[j] + '.png'
            plt.savefig(fname)

        if 0:   #Not doing the plot over a map now
            #raw_input(' hit return for contourf plot over map for this PBAR ')
            plt.clf()
            #extent = (360-124.23, 360-124.16, 41.73,41.785)
            #extent = (235.77, 235.85, 41.7348,41.785)
            plt.imshow(CCmap,extent=region_extent)
            plt.contourf(xp[:-1,:-1],yp[:-1,:-1],zetaplot+Bplot,clines,\
                     colors=colors2,alpha=0.8)
            cbar = plt.colorbar()
            cbar.set_ticks(clines)
            #clines3 = [1e-3,1,2]
            #plt.contour(xp[:-1,:-1],yp[:-1,:-1],zetaplot,clines3,colors='w')
            plt.contour(xp[:-1,:-1],yp[:-1,:-1],B,[0],colors='w')
            plt.ticklabel_format(format='plain',useOffset=False)
            plt.xticks(rotation=20)
            plt.axis(region_extent)
            a = plt.gca()
            a.set_aspect(1./cos(aspectNO*pi/180.))

            titlestr='eta-contours (h+B0) for p= ' + str(pbarlist[j])
            plt.title(titlestr)
            fname=Eventpath + '/eta-contours_prob_' + titlelist[j] + '_map.png'
            plt.savefig(fname)

    outfile.close()
    #raw_input(' hit return to end session ')
    return zetapbar

    
if __name__ == '__main__':
    ###### TYPICAL DATA  #########################################################
    projectdir = os.environ['FEMA']
    scenarioname = 'scenario_test2'
    scenariodir = projectdir + '/' + scenarioname
    fg_types_todo_file = scenariodir + '/fg_types_todo.txt'
    fg_2Dgrids_todo_file = scenariodir + '/fg_2Dgrids_todo.txt'

    print 'INSIDE if name is main, data echo: '
    print ' '
    print 'The project directory was: ',projectdir
    print ' '
    print 'The scenario directory was: ',scenariodir
    print ' '
    
    fg_types_todo = loadtxt(fg_types_todo_file,dtype=integer)
    print 'fg_types_todo where 1 means transects, 2 means 2Dgrids is: '
    print fg_types_todo
    print ' '

    depth=True; speed=False; mflux=False;

    if (depth == True):  #using depth only now
        zeta = linspace(0,12,121)
        whichzetas=[0,5,10,15,20,25,30]  #zeta=0.0,0.5,1.0,1.5,2.0,2.5,3.0
                                         #for p-contours
        whichzetas=array(whichzetas)

    elif (speed==True): #fix later if used
        zeta1=[0.0,.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5]
        zeta2=[10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0]
        zeta3=[14.5,15.0,15.5,16.0,16.5,17.0]
        zeta=array(zeta1+zeta2+zeta3)
        whichzetas=[0,2,4,6,8,10,12,14,16,18,20]  #zeta=0.0,1,2,3,4,5,6,7,8,9,10 m/sec.
        whichzetas=array(whichzetas)

    else: #mflux is true, fix later if used
        #zeta1=[0.0,.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5]
        #zeta2=[10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0]
        #zeta3=[14.5,15.0,15.5,16.0,16.5,17.0]
        #zeta=array(zeta1+zeta2+zeta3)
        #zeta=zeta*100
        zeta=[0.01, 0.05] + list(linspace(0.1,2,20)) + [4,6,8,10,20,50,100,200,400,600,1000,1500,1700]
        zeta=array(zeta)
        whichzetas=[2,11,25,28]
        whichzetas=array(whichzetas)

    print 'The whichzetas for the p-contours were: '
    print whichzetas
    print ' '

    pbarlist=[.1,.2,.3,.4,.5,.6,.7,.8,.9]   #9 zeta-contours
    print 'The pbarlist for the 9 zeta-contours was: '
    print pbarlist
    print ' '

    ######  Now loop over all the work to be done, and call hazard_curves_2Dgrid_plot ###
    ######  for each fixed grid                                                       ###
    for itype in fg_types_todo:
        if (itype == 2):            #Do all the 2D grid type of fixed grids

            #Read the fixed_big_grid_file only once. It is the same little plot for all
            fixed_big_grid_file = projectdir + '/DataFiles/' + 'CCtopo_1_3sec.tt3'
            print 'The fixed grid file for the biggest fixed grid region used for plotting '
            print 'little map in upper right graph corner: '
            print fixed_big_grid_file
            print ' '
            from clawpack.geoclaw.topotools import Topography
            CCtopo = Topography()
            CCtopo.read(fixed_big_grid_file,3)

            fg_2Dgrids_todo=loadtxt(fg_2Dgrids_todo_file,dtype=integer)
            if (size(fg_2Dgrids_todo) == 1):
                fg_2Dgrids_todo = array([fg_2Dgrids_todo])
            for iwhich in fg_2Dgrids_todo:
                outdir = scenariodir + '/' + scenarioname + '_2Dgrid_fg' + str(iwhich)
                hazard_curves_file = outdir + '/hazardcurves_probs.npy'
                fixed_grid_file = projectdir+'/fine_runs/'+'fine_xyB0_fg'+str(iwhich) +'.npy'

                ###Note, need to read Npts1 and Npts2 from some file region_dimensions.txt
                ###      in outdir that make_scenerio.py put there. 
                ###      The file has length two with the values  Npts1 and  Npts2

                dimension_file=outdir + '/' + 'region_dimensions.txt'
                dimensions=loadtxt(dimension_file,dtype=integer)
                Npts1=dimensions[0]; Npts2=dimensions[1];
                print 'The 2Dgrid was Npts1 x Npts2 points in longxlat dimensions'
                print Npts1,Npts2
                print ' '

                #Set this automatically within the plotter for these tests
                #Set whichrows in the plotter
                whichrows=None

                print '#####################'
                print 'CALLING hazard_curves_plot for fixed grid number: ',iwhich
                print 'The fixed grid xyB0 file was: ',fixed_grid_file
                print ' '
                print 'The hazard curves .npy file was: '
                print hazard_curves_file
                print ' '

                fixed_big_grid_file = projectdir + '/DataFiles/' + 'CCtopo_1_3sec.tt3'
                print 'The fixed grid file for the biggest fixed grid region used for plotting '
                print 'little map in upper right graph corner: '
                print fixed_big_grid_file
                print ' '
                from clawpack.geoclaw.topotools import Topography
                CCtopo = Topography()
                CCtopo.read(fixed_big_grid_file,3)

                zetapbar=hazard_curves_2Dgrid_plot(hazard_curves_file,fixed_grid_file,\
                              zeta,whichzetas,outdir,pbarlist,Npts1,Npts2,\
                              CCtopo.X,CCtopo.Y,CCtopo.Z,\
                              whichrows=whichrows,speed=speed,depth=depth,mflux=mflux)
                #Write zetapbar to a file.
                if 0:
                    if (depth==True):
                        outZetafile=outdir + '/' + 'ZetaFloods.npy'
                    elif (speed==True):
                        outZetafile=outdir + '/' + 'SpeedFloods.npy'
                    else:
                        outZetafile=outdir + '/' + 'MfluxFloods.npy'
                    save(outZetafile,zetapbar)
                print ' '
