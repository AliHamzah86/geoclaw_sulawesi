"""

Usage: 
 >>>
  from hazard_curves_transect_plot import hazard_curves_transect_plot

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
                   columns in zeta for which we want transect
                   plots of probabilities. For example, 
                   whichzetas=array([0,3]) would use the
                   probabilities corresponding to P1(col 0),
                   and P4(col 3).  Each of the whichzetas is
                   an exceedance value (of depth). So
                   plotting what is in a column is a
                   p-transect, for example, a representation 
                   of the probabilities at all the points of
                   exceeding a whichzeta value (say 2 meters).
                   

      Eventpath  = directory path where to save the plots

      pbarlist: 
             Made this an input, rather than set within the code.  It is
             a list of "floods" for the zeta-transects.
             For example, pbarlist=[.01,.002] would be the 100 year flood
             or the 500 year flood.  Another example, if we are only interested
             in conditional probabilities (not whether or not Cascadia happens,
             for example), we might take pbarlist=[.3,.5,.7,.8,.9].  The
             entries are the probability that a level is exceeded. 

      X_CC,Y_CC,Z_CC:
             These are 2D structures for the little contour plot in upper right
             corner of some plots.

      Example:
         import os
         CCdir = os.environ['CCdir']
         probsnpy='hazardcurves_probs.npy'
         Testsdir= CCdir + '/Realizations_Test/'
         fixed_grid_file = CCdir + '/Datafiles/fixedgrid_xyB_test.npy'

         hazard_curves_file  = Testsdir +'EVENTeall' + '/'+ probsnpy 

         See the __main__ below for setting zeta

         whichrows is now set within the routine below based on the number
                   of grid points on the transect - done automatically.
                   Particular gauges could be added in the calling program
                   by setting whichrows as a list with these gauge indices,
                   then this program appends to the inputted list, the
                   automatic numbers.

         speed=True if we are doing plots of speeds; otherwise speed=False.
         depth=True if we are doing plots of depth; otherwise depth=False.
         mflux=True if we are doing plots of mflux; otherwise mflux=False.

    OUTPUT:  Does a plot for each row of the hazard_curves_file that
             is indicated in whichrows.  The first row is row 0.
             User may be prompted to step through the plots.

             Also does zeta-transect plots where p is specified.
             Say the p=.01 "event" (100 year flood).
             Also could write the Zeta values for the p's specified.

             For example:

             outZetafile=Testsdir + 'EVENTeall' + '/' + 'Zeta100_500.npy'

             Also p-transect plots of the probabilities
             associated with a particular zetai level as specified in
             the list whichzetas. User may be prompted to step through the
             plots. The list whichzetas are the columns in the 
             hazard_curves_file that we want: 

             For example:

             whichzetas=[0,5]

             would correspond to zeta0 and zeta5. 


"""
from numpy import array,shape,load,nan,nanmax,isnan,linspace,sum,where
from numpy import zeros,ones,ma,amin,amax,logical_or,logical_and
from numpy import int,integer,ceil
from numpy import loadtxt,size
from numpy import cos,pi,arange
from matplotlib.mlab import find
from matplotlib import colors
import matplotlib.pyplot as plt
import os,sys

def hazard_curves_transect_plot(hazard_curves_file,fixed_grid_file,\
       zeta,whichzetas,Eventpath,pbarlist,\
       X_CC,Y_CC,Z_CC,\
       whichrows=None,speed=False,depth=True,mflux=False):

    aspectNO=41.75                                     #aspect latitude for plotting
                                                       #specific to Crescent City
    commentsfile=Eventpath + '/hc_plot_run.output'
    outfile=open(commentsfile,'w')

    fprobs=load(hazard_curves_file)
    d=load(fixed_grid_file)
    x=d[:,0]
    y=d[:,1]
    B=d[:,2]           #This is B0 for the FEMA project

    #Plot transect on the contour map
    plt.figure(1)
    plt.clf()
    plt.contour(X_CC, Y_CC, Z_CC, arange(0,21,2), colors='g')
    plt.xticks(rotation=20)
    plt.xlim(-124.215,-124.18)
    plt.ylim(41.735,41.768)
    plt.ticklabel_format(format='plain',useOffset=False)
    a = plt.gca()
    a.set_aspect(1./cos(aspectNO*pi/180.))
    plt.hold(True)
    plt.plot(x,y,'r-')
    plt.title('Transect on contour map')
    fname=Eventpath + '/Transect_on_contour_map.png'
    plt.savefig(fname)
    plt.close(1)

    #Determine where the first land points (includes B0=0) start on the transect 
    land_indices = find(B >= 0.0)              #includes B0=0 points
    B_land=ma.masked_where(B < 0.0,B)         #includes B0=0 points, B_land same size as B
    upper2=nanmax(B_land)

    y1=y[land_indices[0]]; y2=y[land_indices[-1]];
    print ' '
    print 'The first and last land points (B0>=0) are at latitudes: '
    print y1,y2
    print ' '


    Npts=len(x)
    outfile.write(' THE number of points on this transect fixed grid is : ')
    outfile.write('%s ' %Npts)
    outfile.write(' \n')
    outfile.write(' \n')

    if (whichrows == None):
        whichrows=[]

    jj=int(Npts/5.0)
    for j in [1,2,3,4,5]:
        whichrowsindex=j*jj -1    #python starts at 0
        whichrows.append(whichrowsindex)
    whichrows=array(whichrows)

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

    outfile.write(' The location indices for hazard plots are: \n')
    outfile.write('%s ' %whichrows)
    outfile.write(' \n')
    outfile.write(' \n')

    outfile.write(' The indices for zeta-exceedance levels for p-transects are: \n')
    outfile.write('%s ' %whichzetas)
    outfile.write(' \n')
    outfile.write(' \n')

    outfile.write(' The flood levels (p values) for zeta-transects are: \n')
    outfile.write('%s ' %pbarlist)
    outfile.write(' \n')
    outfile.write(' \n')

    plt.figure(11,(12,8))
    fprobs121=zeros(121)        
    for irow in whichrows:       #Grid Point Hazard Curve Plots
        plt.clf()
        fprobs121=fprobs[irow,:]
        masking=where(isnan(fprobs121),True,False)
        fprobs121=ma.masked_where(masking,fprobs121)

        #fprobs121=where(fprobs121 <= 1e-10,1e-10,fprobs121) #So semilog can work

        if (mflux == True):
            fprobs121=where(fprobs121 <= 1e-10,1e-10,fprobs121) #mflux use loglog
            plt.loglog(zeta,fprobs121,'ro')
            plt.hold(True) 
            plt.loglog(zeta,fprobs121,'r-')
        else:                                    #Not doing semilog for depth or speed plots 
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
        plt.title('At Longitude, Latitude, Bathy:  (%9.4f,%9.4f,%9.4f) ' \
                                                      %(x[irow],y[irow],d[irow,2]))

        ######### putting little map on the plot, showing the location
        #fixed grid was a 1D transect
        plt.axes([.7,.65,.25,.25])  #.7,.65 is lower left ratio of figure,.25,.25 width,height
        plt.contour(X_CC, Y_CC, Z_CC, arange(0,21,2), colors='k')
        a = plt.gca()
        a.set_aspect(1./cos(aspectNO*pi/180.))

        plt.plot(x,y,'r-')
        xi=x[irow]; yi=y[irow];
        plt.plot([xi],[yi],'ro')
        plt.xticks([])
        plt.yticks([])
        #########

        #raw_input(' hit return to get the next plot ')
        plt.savefig(fname)
    plt.close(11)

    plt.figure(12,(12,8))         #one of these for all whichzeta
    plt.clf()      
    plt.figure(15,(12,8))         #one of these for each whichzeta
    plt.clf()

    ymid=(y2+y1)/2.0; y1str=str(y1); y2str=str(y2); ymstr=str(ymid);
    ystart_str=str(y[0])
    xeps=-0.001

    legendlist=[]
    for izeta in whichzetas:      #p-transect plots for heights in whichzetas
        legendlist.append('zeta = %5.2f m.' %zeta[izeta])
        zetarest=fprobs[:,izeta]
        #masking=logical_or(isnan(zetarest),zetarest <=0.0)
        masking=isnan(zetarest)
        zetarest=ma.masked_where(masking,zetarest)

        #raw_input(' hit return for p-transect for zeta= %5.2f ' %(zeta[izeta]))
        plt.figure(15)
        plt.plot(y,zetarest,'r-')
        #plt.xlim(y1+xeps,y2); plt.ylim(-0.1,1.5);
        plt.xlim(y[0],y[-1]); plt.ylim(-0.1,1.1);
        #plt.yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
        plt.yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],\
                   ['0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7',\
                    '0.8','0.9','1.0'])
        plt.xticks([y[0],y1,ymid,y2],[ystart_str,y1str,ymstr,y2str],rotation=20)
        #plt.xticks([y1,ymid,y2],[y1str,ymstr,y2str],rotation=20)
        plt.xlabel('Latitude on the Transect')    #if transect is horizontal, in trouble!
        plt.ylabel('Probability of Exceedance')

        nn=str(int(zeta[izeta]*100)).zfill(5)
        if (depth==True):
            plt.title('p-transect for zeta=%5.2f m' %(zeta[izeta]), fontsize=20)
            fname=Eventpath + '/p-transect_zeta_%s.png' %nn
        elif (speed==True):
            plt.title('p-transect for speed=%5.2f m/sec' %(zeta[izeta]), fontsize=20)
            fname=Eventpath + '/p-transect_speed_%s.png' %nn
        else:
            plt.title('p-transect for momentum flux=%5.2f ' %(zeta[izeta]), fontsize=20)
            fname=Eventpath + '/p-transect_mflux_%s.png' %nn
        plt.savefig(fname)
        plt.clf()           ##clears figure 15

        plt.figure(12)
        plt.plot(y,zetarest)

    plt.close(15)

    #plt.xlim(y1+xeps,y2); plt.ylim(-0.1,1.1);
    plt.xlim(y[0],y[-1]); plt.ylim(-0.1,1.8);
    plt.yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],\
               ['0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0'])
    #plt.xticks([y1,ymid,y2],[y1str,ymstr,y2str],rotation=20)
    plt.xticks([y[0],y1,ymid,y2],[ystart_str,y1str,ymstr,y2str],rotation=20)
    plt.legend(legendlist,loc='upper right')
    plt.xlabel('Latitude on the Transect')    #if transect is horizontal, in trouble!
    plt.ylabel('Probability of Exceedance')
    plt.title('p-transects for Depth ', fontsize=20)
    fname=Eventpath + '/p-transect_all_zetas.png'
    plt.savefig(fname)
    plt.close(12)

    #Determine the zetapt01 and zetapt002 for the 100 yr and 500 yr events 
    #and store in the two columns of zetapbar
    #Find zetapt01 and zetapt002, zeta300, zeta975, zeta2500, zeta5000

    #Note, made the program more general, the particular "floods" are in pbarlist below
    #and are not limited to the yr events stated above.

    sh=shape(fprobs)
    Npts=sh[0]
    Nzeta=sh[1]

    numpbar = len(pbarlist)
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

        outfile.write('The number of grid points was: %5d \n' %Npts)
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

    #Note: okval is a list of True or False for whether there is
    #      a hazard flood for this value or not.  
    #      So need to check the okval list to know.
    #

    uppersave=[]; legendlist=[];
    B0_neg = where(B < 0.0, B, nan)
    B0_pos = where(B < 0.0, 0, B)          #used to make eta. If water point, B already added to h
                                           #for land points, zetaplot only has h effect, not B.
                                           #We will add B0_pos to zetaplot

    plt.figure(12,(12,8))                  #all on same plot, no legend
    plt.clf()
    plt.fill_between(y,B0_neg,0,color=[.6,.6,1])
    plt.plot(y,B,'g-',linewidth=3)                   

    plt.figure(13,(12,8))                  #all on same plot, with legend
    plt.clf()
    legendlist.append('topo(B0)')
    plt.fill_between(y,B0_neg,0,color=[.6,.6,1])
    plt.plot(y,B,'g-',linewidth=3)                   

    plt.figure(14,(12,8))                  #each pbar individually
    plt.clf()
    for j in range(len(pbarlist)):
        outfile.write(' \n')
        outfile.write(' PBAR =: %s ' %pbarlist[j])
        zetaplot=zetapbar[:,j]
        #masking=logical_or(isnan(zetaplot),zetaplot <=0.0)
        masking=isnan(zetaplot)
        zetaplot=ma.masked_where(masking,zetaplot)
        Bplot=ma.masked_where(masking,B0_pos)    #mask Bplot at same place as zetaplot
                                                 #Bplot can be masked differently each
                                                 #pass thru the loop, so call it Bplot
                                                 #For example, zetaplot and Bplot are masked
                                                 #where no flood for this probability
                                                 #zeta=0.0 was not even EXCEEDED.
                                                 #Bplot will be 0 for a water point as nothing
                                                 #needs to be added to zetaplot, since zetaplot
                                                 #there already represents B0+h.

        if (okval[j] == True):
            #outfile.write('\n')
            #outfile.write(' FOR PBAR == : %s ' %pbarlist[j])
            outfile.write('max value in zetaplot %7.3f \n' \
                      %(amax(zetaplot)))
            outfile.write('min value in zetaplot %7.3f \n' \
                      %(amin(zetaplot)))
            outfile.write('\n')
            outfile.write('\n')

        else:
            outfile.write('\n')
            outfile.write(' FOR PBAR == : %s ' %pbarlist[j])
            outfile.write('No zeta-transect exists for this PBAR  \n')
            outfile.write('\n')
            outfile.write('\n')

        if 1:
            #raw_input(' hit return for this PBAR transect plot ')
            plt.figure(14)
            plt.subplot(2,1,1)        #plot zeta alone on the first subplot
            #plt.fill_between(y,B0_neg,0,color=[.6,.6,1])
            #plt.plot(y,B,'g-',linewidth=3)
            plt.plot(y,zetaplot,'r-')
            #plt.plot(y,B_land,'g-',linewidth=3)
            #plt.xlim(y1+xeps,y2)

            plt.xlim(y[0],y[-1])
            upperlim=max(18.0,upper2)
            upp=ceil(upperlim)
            uppint=int(upp)
            plt.ylim(-10.0,upp)
            if (upperlim <= 18.0):
                plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0],\
                           ['-10','-5','0','2','4','6','8','10','12','14','16','18'])
            else:
                plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,upp],\
                           ['-10','-5','0','2','4','6','8','10','12','14','16','18',uppint])
            plt.xticks([])
            #plt.xticks([y1,ymid,y2],[y1str,ymstr,y2str],rotation=20)
            #plt.xlabel('Latitude on the Transect')  #Note, in trouble if horizontal transect
            plt.ylabel('meters')

            if (depth==True):
                titlestr='zeta-transect for p= ' + str(pbarlist[j]) + ' with topo(B0)'
                plt.title(titlestr)
                fname=Eventpath + '/zeta-transect_prob_' + titlelist[j] + '.png'
            elif (speed==True):
                titlestr='speed-transect for p= ' + str(pbarlist[j])
                plt.title(titlestr)
                fname=Eventpath + '/speed-transect_prob_' + titlelist[j] + '.png'
            else:
                titlestr='mflux-transect for p= ' + str(pbarlist[j])
                plt.title(titlestr)
                fname=Eventpath + '/mflux-transect_prob_' + titlelist[j] + '.png'

            #plt.savefig(fname)   #saves along with the next subplot
            #plt.clf()            #clears figure 14, now making eta a subplot

        todo=1                    #do the eta plot
        if (logical_and(depth==True,todo==1)):
            #Now put the second subplot on this figure 14, which is the eta, B, and fill
            plt.subplot(2,1,2)
            plt.fill_between(y,B0_neg,0,color=[.6,.6,1])
            plt.plot(y,B,'g-',linewidth=3)
            #plt.plot(y,zetaplot,'r-')       #don't plot zeta on this plot
            plt.plot(y,zetaplot+Bplot,'b-')
            #plt.plot(y,B_land,'g-')  
            #plt.xlim(y1+xeps,y2)

            plt.xlim(y[0],y[-1])

            #We want to continue plotting B0 beyond where zetaplot+Bplot is masked.
            upper1=nanmax(zetaplot+Bplot)
            print ' '
            print ' PBAR =:  ', pbarlist[j]
            print 'The max of h+B0,B0 on transect were: ',upper1,upper2
            print ' '
            upper=max(upper1,upper2)
            uppersave.append(upper)

            upp=ceil(upper)
            plt.ylim(-10.0,upp)
            if (upper <= 18.0):
                plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0],\
                           ['-10','-5','0','2','4','6','8','10','12','14','16','18'])
            else:
                lasttick=str(int(upp))
                plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,upp],\
                           ['-10','-5','0','2','4','6','8',10,'12','14','16','18',lasttick])
            #plt.xticks([y1,ymid,y2],[y1str,ymstr,y2str],rotation=20)
            plt.xticks([y[0],y1,ymid,y2],[ystart_str,y1str,ymstr,y2str],rotation=20)
            plt.xlabel('Latitude on the Transect')  #Note, in trouble if horizontal transect
            plt.ylabel('Blue: Surface(h+B0), Green: Bathymetry(B0)')

            titlestr='eta-transect for p= ' + str(pbarlist[j]) + ' with topo(B0)'
            plt.title(titlestr)
            fname=Eventpath + '/zeta-eta-transects_prob_' + titlelist[j] + '.png'
            plt.savefig(fname)
            plt.clf()          ##clears figure 14

            legendlist.append('p = %5.2f' %pbarlist[j])
            plt.figure(13)
            plt.plot(y,zetaplot+Bplot)
            plt.figure(12)
            plt.plot(y,zetaplot+Bplot)
    
    plt.close(14)
    if (logical_and(depth==True,todo==1)):  #finish figure(12) and figure(13) of all eta plot
        uppersave=array(uppersave)
        maxupper=max(uppersave)
        upp=ceil(maxupper)
        lasttick=str(int(upp))

        plt.figure(12)
        plt.xlim(y[0],y[-1])
        plt.ylim(-10.0,upp)
        if (maxupper <= 18.0):
            plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0],\
                       ['-10','-5','0','2','4','6','8','10','12','14','16','18'])
        else:
            plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,upp],\
                       ['-10','-5','0','2','4','6','8',10,'12','14','16','18',lasttick])
        #plt.xticks([y1,ymid,y2],[y1str,ymstr,y2str],rotation=20)
        plt.xticks([y[0],y1,ymid,y2],[ystart_str,y1str,ymstr,y2str],rotation=20)
        #plt.legend(legendlist,loc='lower right')
        plt.xlabel('Latitude on the Transect')  #Note, in trouble if horizontal transect
        plt.ylabel('meters')
        titlestr='Surface(h+B0) elevation along transect with topo(B0)'
        plt.title(titlestr)
        fname=Eventpath + '/eta-transects_all_nolegend' + '.png'
        plt.savefig(fname)
        plt.close(12) 

        plt.figure(13)
        plt.xlim(y[0],y[-1])
        plt.ylim(-10.0,upp)
        if (maxupper <= 18.0):
            plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0],\
                       ['-10','-5','0','2','4','6','8','10','12','14','16','18'])
        else:
            plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,upp],\
                       ['-10','-5','0','2','4','6','8',10,'12','14','16','18',lasttick])
        #plt.xticks([y1,ymid,y2],[y1str,ymstr,y2str],rotation=20)
        plt.xticks([y[0],y1,ymid,y2],[ystart_str,y1str,ymstr,y2str],rotation=20)
        plt.legend(legendlist,loc='lower right')
        plt.xlabel('Latitude on the Transect')  #Note, in trouble if horizontal transect
        plt.ylabel('meters')
        titlestr='Surface(h+B0) elevation along transect with topo(B0)'
        plt.title(titlestr)
        fname=Eventpath + '/eta-transects_all' + '.png'
        plt.savefig(fname)
        plt.close(13) 
        
        
    outfile.close()
    #raw_input(' hit return to end session ')
    return zetapbar

    
if __name__ == '__main__':
    ###### TYPICAL DATA  #########################################################

    projectdir = os.environ['FEMA']
    scenarioname = 'scenario_test2'
    scenariodir = projectdir + '/' + scenarioname
    fg_types_todo_file = scenariodir + '/fg_types_todo.txt'
    fg_transects_todo_file = scenariodir + '/fg_transects_todo.txt'
    fixed_big_grid_file = projectdir + '/DataFiles/' + 'CCtopo_1_3sec.tt3'

    print 'INSIDE if name is main, data echo: '
    print ' '
    print 'The project directory was: ',projectdir
    print ' '
    print 'The scenario directory was: ',scenariodir
    print ' '
    print 'The whole region grid xyB0 file used for little contour plots: '
    print 'in the upper right corner of some plots was: '
    print fixed_big_grid_file
    print ' '
    
    fg_types_todo = loadtxt(fg_types_todo_file,dtype=integer)
    print 'fg_types_todo where 1 means transects, 2 means 2Dgrids is: '
    print fg_types_todo
    print ' '

    depth=True; speed=False; mflux=False;

    if (depth == True):  #using depth only now
        zeta = linspace(0,12,121)
        whichzetas=[0,10,20,30]  #zeta=0.0,1.0,2.0,3.0
                                 #for p-transects
        whichzetas=array(whichzetas)

    elif (speed==True): #fix later if used
        zeta1=[0.0,.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5]
        zeta2=[10.0,10.5,11.0,11.5,12.0,12.5,13.0,13.5,14.0]
        zeta3=[14.5,15.0,15.5,16.0,16.5,17.0]
        zeta=array(zeta1+zeta2+zeta3)
        whichzetas=[0,2,4,6,8,10,12,14,16,18,20]  #zeta=0.0,1,2,3,4,5,6,7,8,9,10 m/sec.
        whichzetas=array(whichzetas)

    else: #mflux is true, fix later if used
        zeta=[0.01, 0.05] + list(linspace(0.1,2,20)) + [4,6,8,10,20,50,100,200,400,600,1000,1500,1700]
        zeta=array(zeta)
        whichzetas=[2,11,25,28]
        whichzetas=array(whichzetas)

    print 'The whichzetas for the p-transects were: '
    print whichzetas
    print ' '

    pbarlist=[.1,.2,.3,.4,.5,.6,.7,.8,.9]   #9 zeta-transects
    print 'The pbarlist for the 9 zeta-transects was: '
    print pbarlist
    print ' '

    ######  Now loop over all the work to be done, and call hazard_curves_transect_plot ###
    for itype in fg_types_todo:
        if (itype == 1):            #Do all the transect type of fixed grids

            # Read the following once as plots use the same little contour
            print 'The whole region grid xyB0 file used for little contour plots: '
            print 'in the upper right corner of some plots was: '
            print fixed_big_grid_file
            print ' '
            from clawpack.geoclaw.topotools import Topography
            CCtopo = Topography()
            CCtopo.read(fixed_big_grid_file,3)

            fg_transects_todo=loadtxt(fg_transects_todo_file,dtype=integer)
            if (size(fg_transects_todo) == 1):
                fg_transects_todo = array([fg_transects_todo])
            for iwhich in fg_transects_todo:
                outdir = scenariodir + '/' + scenarioname + '_transect_fg' + str(iwhich)
                hazard_curves_file = outdir + '/hazardcurves_probs.npy'
                fixed_grid_file = projectdir+'/fine_runs/'+'fine_xyB0_fg'+str(iwhich) +'.npy'

                #Set this automatically within the plotter for this test
                whichrows=None     

                print '#####################'
                print 'CALLING hazard_curves_transect_plot with fixed grid number: ',iwhich
                print 'The fixed grid xyB0 file was: ',fixed_grid_file
                print ' '
                print 'The hazard curves .npy file was: '
                print hazard_curves_file
                print ' '

                zetapbar=hazard_curves_transect_plot(hazard_curves_file,fixed_grid_file,\
                          zeta,whichzetas,outdir,pbarlist,\
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
            
