"""

Usage: 
 >>>
  from hazard_curves_transect_compare import hazard_curves_transect_compare

    where for example

    INPUT: 
      hazard_curves_file = the .npy file used for the plots:
         It is a directory of .npy files used for comparisons of the
         scenarios on this particular transect and has the format:

           P1 P2 P3 PNzeta 
        
      OutZetafile = the .npy file of the zetapbar values produced by hc_plot_run.py
         It is a directory of .npy files used for the "flood" comparsions where the flood
         level is given by the whichzeta list.

      fixed_grid_file = file with x y B in colms, a .npy file
        (for the FEMA project, B is B0, initial bathy)
        Will change in the calling program for each transect.
        Fixed in this program. Made this a list, so the things
        being compared could be say a fine resolution grid and a
        coarse resolution grid.  This would require two different B
        values with the same x and y.  To know which scenario has
        which resolution, we added the scenarioresolution list below.

      scenarioresolution = ['fine','fine'] to compare two fine scenarios, or
                         = ['fine','coarse'] if the two are fine and coarse.

      zeta = numpy array of Nzeta exceedance values 

      whichrows = numpy array of integers representing the
                  rows of the hazard_curves_file to plot
                  the first row in numpy array is index 0!
                  Each row represents one hazard curve at a
                  particular point. Will be determined automatically
                  for this transect in this program.
                  It is now set within the routine below based on the number
                  of grid points on the transect - done automatically.
                  Particular gauges could be added in the calling program
                  by setting whichrows as a list with these gauge indices,
                  then this program appends to the inputted list, the
                  automatic numbers.

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
                   
      Eventpath  = directory path where to save the comparison plots 
                   for this particular transect. 

      legendnames= Names to use on the comparison plots for the scenarios

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
         Testsdir1= CCdir + '/Realizations_Test1/'
         Testsdir2= CCdir + '/Realizations_Test2/'
         fixed_grid_file = [CCdir + '/Datafiles/fixedgrid_xyB_test.npy',
                            CCdir + '/Datafiles/fixedgrid_xyB_test.npy']

         hazard_curves_file  = [Testsdir1 +'EVENTeall' + '/'+ probsnpy, 
                               [Testsdir2 +'EVENTeall' + '/'+ probsnpy] 

         OutZetafile = [Testsdir1 + 'EVENTeall' + '/' + 'Zeta100_500.npy',
                        Testsdir2 + 'EVENTeall' + '/' + 'Zeta100_500.npy]

         Eventpath = CCdir + 'Realizations/comparisons/comparisons_1/transect_fg1'

         See the __main__ below for setting zeta

         speed=True if we are doing plots of speeds; otherwise speed=False.
         depth=True if we are doing plots of depth; otherwise depth=False.
         mflux=True if we are doing plots of mflux; otherwise mflux=False.

    OUTPUT:  Hazard Curves comparison plots for the whichrows grid points.
             p-transects comparisons plots for various zeta levels in whichzetas.
             zeta-transects comparisons plots for various p levels in pbarlist. 

             These comparison plots are all stored in the directory Eventpath
             (which is called comparisondir in the calling program). A commentsfile
             of the job run is also stored there.


"""
from numpy import array,shape,load,nan,nanmax,isnan,linspace,sum,where
from numpy import zeros,ones,ma,amin,amax,logical_or,logical_and
from numpy import int,integer,ceil,floor
from numpy import loadtxt,size
from numpy import cos,pi,arange
from matplotlib.mlab import find
from matplotlib import colors
import matplotlib.pyplot as plt
import os,sys

def hazard_curves_transect_compare(hazard_curves_file,fixed_grid_file,\
       Outzetafile,\
       zeta,whichzetas,Eventpath,pbarlist,\
       X_CC,Y_CC,Z_CC,legendnames,scenarioresolution,\
       whichrows=None,speed=False,depth=True,mflux=False):

    aspectNO=41.75                                     #aspect latitude for plotting
                                                       #specific to Crescent City
    commentsfile=Eventpath + '/hc_compare_run.output'
    outfile=open(commentsfile,'w')

    fprobs_list=[]; zetapbar_list=[];
    upper2_list=[]; y1_list=[]; 
    B0_neg_list=[]; B0_pos_list=[]; B_list=[];
    fine_found=0; use_index=0;                    #Use the first scenario for the bathy plot when
                                                  #no fine bathy is found below.

    plt.figure(10,(12,8))                         #Plot all the bathy from the scenarios
    plt.clf()
    num_scenarios=len(hazard_curves_file)
    minB=0.0
    for j in range(num_scenarios):
        d=load(fixed_grid_file[j])
        x=d[:,0]; y=d[:,1];                       #long and lat for this scenario
        B=d[:,2]                                  #This is B0 for this scenario for the FEMA project

        if (min(B) < minB):
            minB=min(B)

        plt.plot(y,B,'-')                         #Plot this scenario bathymetry
        B_list.append(B)
        B0_neg_list.append(where(B<0.0, B,nan))   #Used for plotting
        B0_pos_list.append(where(B < 0.0, 0, B))  #Used to make eta. If water point, B already added
                                                  #to h.
        land_indices = find(B >=0.0)
        B_land=ma.masked_where(B < 0.0,B)         #includes B0=0 points, B_land same size as B
        upper2_list.append(nanmax(B_land))             #find plotting limits for bathymetry
        y1_list.append(y[land_indices[0]]);       #first land point of this scenario
        outfile.write(' \n')
        outfile.write(' For Scenario: %s ' %legendnames[j])
        outfile.write(' \n')
        outfile.write(' The first land point (B0>=0) is at latitude:  \n')
        outfile.write(' %s ' %y1_list[j])
        outfile.write(' \n')
        outfile.write(' \n')
        fprobs_list.append(load(hazard_curves_file[j]))
        zetapbar_list.append(load(Outzetafile[j]))
        if (logical_and(scenarioresolution[j] == 'fine',fine_found==0)):
                    fine_found=1; use_index=j;    #using the first fine one found

    #All the scenarios for this transect have the same x and y,
    #so x and y are now set to those of the last scenario which is fine
  
    upper2=max(upper2_list)
    outfile.write(' upper Bathy on land over all scenarios was: %s ' %upper2)
    outfile.write(' \n')

    #Make ticks at ystart, y1, ymid, y2
    lat_eps=0.001/2.0
    y1=min(y1_list); y2=y[-1];
    if (abs(y1 - max(y1_list)) <= lat_eps):
        ymid=(y2+y1)/2.0
    else:
        ymid=max(y1_list)
    y1str=str(y1); y2str=str(y[-1]); ymstr=str(ymid);
    ystart_str=str(y[0])

    ##Finish the bathymetry-transect comparison plot
    plt.xlim(y[0],y[-1])
    upp=ceil(upper2)
    plt.ylim(-10.0,upp)
    if (upper2 <= 18.0):
        firsttick=str(int(floor(minB)))
        plt.yticks([floor(minB),-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0],\
                   [firsttick,'-5','0','2','4','6','8','10','12','14','16','18'])
    else:
        firsttick=str(int(floor(minB)))
        lasttick=str(int(upp))
        plt.yticks([floor(minB),-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,'18',upp],\
                   [firsttick,'-5','0','2','4','6','8',10,'12','14','16','18',lasttick])
    plt.xticks([y[0],y1,ymid,y2],[ystart_str,y1str,ymstr,y2str],rotation=20)
    plt.xlabel('Latitude on the Transect')  #Note, in trouble if horizontal transect
    plt.ylabel('Bathymetry in meters')
    titlestr='Bathymetries for all Scenarios '
    plt.title(titlestr)
    plt.legend(legendnames,loc='upper left')
    fname=Eventpath + '/bathymetry-transects.png'
    plt.savefig(fname)
    plt.close(10)


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
        #plt.title('At Longitude, Latitude, Bathy:  (%9.4f,%9.4f,%9.4f) ' \
        #                                              %(x[irow],y[irow],d[irow,2]))
        plt.title('At Longitude, Latitude (%9.4f,%9.4f) \n' % (x[irow],y[irow]) \
                  +'Scenario Bathymetry Values: ' + num_scenarios*'%9.4f '\
                   %tuple(B_vals))

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

        plt.savefig(fname)    #save plot for this irow in whichrow
        #raw_input(' hit return to get the next plot ')
    plt.close(11)

    ### Now do the p-transect Comparison Plots
    plt.figure(15,(12,8))         #one of these for each whichzeta
    for izeta in whichzetas:      #p-transect plots for heights in whichzetas
        plt.clf()                 #clears figure 15 before next izeta in whichzetas 
        for j in range(num_scenarios):
            fprobs=fprobs_list[j]
            zetarest=fprobs[:,izeta]
            #masking=logical_or(isnan(zetarest),zetarest <=0.0)
            masking=isnan(zetarest)
            zetarest=ma.masked_where(masking,zetarest)

            #raw_input(' hit return for p-transect for zeta= %5.2f ' %(zeta[izeta]))
            plt.plot(y,zetarest)
            plt.hold(True)
        plt.legend(legendnames,loc='upper right')
        #plt.xlim(y1+xeps,y2); plt.ylim(-0.1,1.5);
        plt.xlim(y[0],y[-1]); plt.ylim(-0.1,1.5); 
        #plt.yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
        plt.yticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],\
                   ['0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7',\
                    '0.8','0.9','1.0'])
        plt.xticks([y[0],y1,ymid,y2],[ystart_str,y1str,ymstr,y2str],rotation=20)
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
        plt.savefig(fname)  ##save for this izeta in whichzetas
    plt.close(15)


    ### Now do the "Floods" Comparison Plots
    titlelist=[]
    for j in range(len(pbarlist)):
        titlelist.append('p=' + str(pbarlist[j]))
    
    plt.figure(14,(12,8))                  #each pbar individually
    for j in range(len(pbarlist)):
        outfile.write(' \n')
        outfile.write(' PBAR =: %s ' %pbarlist[j])
        outfile.write(' \n')
        plt.clf()                          #clears figure 14 for next pbar

        plt.subplot(2,1,1)                 #plot zeta alone on first subplot
        for jsc in range(num_scenarios):
            zetapbar=zetapbar_list[jsc]
            zetaplot=zetapbar[:,j]
            #masking=logical_or(isnan(zetaplot),zetaplot <=0.0)
            masking=isnan(zetaplot)
            zetaplot=ma.masked_where(masking,zetaplot)

            #make a zetaplot for this jsc on figure 14
            #raw_input(' hit return for this PBAR transect plot ')
            plt.plot(y,zetaplot,label=legendnames[jsc])
            plt.hold(True)
        #plt.xlim(y1+xeps,y2)
        plt.xlim(y[0],y[-1])
        upperlim=max(18.0,upper2)
        upp=ceil(upperlim)
        plt.ylim(-10.0,upp)
        uppint=str(int(upp))

        #ADDED negative values to include some harbor bathy on plot.
        if (upperlim <= 18.0):
            plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0],\
                       ['-10','-5','0','2','4','6','8','10','12','14','16','18'])
        else:
            plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,upp],\
                       ['-10','-5','0','2','4','6','8','10','12','14','16','18',uppint])
        plt.xticks([])
        plt.ylabel('meters')

        if (depth==True):
            titlestr='zeta-transect for p= ' + str(pbarlist[j])
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
        plt.legend(loc='upper right')
        #plt.savefig(fname)       #not saving zeta plot separately, making a subplot

        #plt.clf()                #clears figure 14, now making eta a subplot

        todo=1                    #do the eta plot on figure 14, as the second subplot
        if (logical_and(depth==True,todo==1)):
            #Now put the second subplot on this figure 14, which is the eta, B, and fill
            #The bathymetry plotted on this plot will be from the fine grid if one of the
            #scenarios was fine.  This scenarios eta value however, will be calculated based
            #on the bathymetry used for each particular scenario.

            plt.subplot(2,1,2)
            plt.fill_between(y,B0_neg_list[use_index],0,color=[.6,.6,1])
            plt.plot(y,B_list[use_index],'g-',linewidth=3)
            uppersave=[]
            for jsc in range(num_scenarios):
                zetapbar=zetapbar_list[jsc]
                B0_pos=B0_pos_list[jsc]                      #B0_pos for this scenario
                zetaplot=zetapbar[:,j]
                #masking=logical_or(isnan(zetaplot),zetaplot <=0.0)
                masking=isnan(zetaplot)
                zetaplot=ma.masked_where(masking,zetaplot)
                Bplot=ma.masked_where(masking,B0_pos) #mask B at same place as zetaplot
                                                      #B can be masked differently each
                                                      #pass thru the loop, so call it Bplot
                                                      #For example, zetaplot and Bplot are masked
                                                      #where no flood for this probability
                                                      #zeta=0.0 was not even EXCEEDED.
                                                      #Bplot will be 0 for a water point as nothing
                                                      #needs to be added to zetaplot, since zetaplot
                                                      #there already represents B0+h.

                #We want to continue plotting B0 beyond where zetaplot+Bplot is masked.
                upper1=nanmax(zetaplot+Bplot)

                outfile.write(' %s: max(h+B0), max(B0) are: ' %legendnames[jsc])
                outfile.write(' %s %s ' %(upper1,upper2_list[jsc]))
                outfile.write(' \n')
                
                upper=max(upper1,upper2)
                upper_use=max(upper,18)
                uppersave.append(upper_use)
                plt.plot(y,zetaplot+Bplot,label=legendnames[jsc])
                plt.hold(True)

            uppersavemax=max(uppersave)
            #plt.xlim(y1+xeps,y2)
            plt.xlim(y[0],y[-1]) 
            upp=ceil(uppersavemax)
            plt.ylim(-10.0,upp)
            #plt.ylim(-0.5,uppersavemax)         #upper limit over all scenarios on the plot
            if (uppersavemax <= 18.0):
                plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0],\
                           ['-10','-5','0','2','4','6','8','10','12','14','16','18'])
            else:
                lasttick=str(int(upp))
                plt.yticks([-10.0,-5.0,0.0,2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0,upp],\
                           ['-10','-5','0','2','4','6','8','10','12','14','16','18',lasttick])
            plt.xticks([y[0],y1,ymid,y2],[ystart_str,y1str,ymstr,y2str],rotation=20)
            plt.xlabel('Latitude on the Transect')  #Note, in trouble if horizontal transect
            plt.ylabel('meters')

            titlestr='eta-transect for p= ' + str(pbarlist[j]) + ' with topo B0'
            plt.title(titlestr)
            plt.legend(loc='lower right')
            fname=Eventpath + '/zeta-eta-transects_prob_' + titlelist[j] + '.png'
            plt.savefig(fname)
    plt.close(14)
        
    outfile.close()
    #raw_input(' hit return to end session ')
    
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
    print 'The scenarionames (also used for plot legend names) was: '
    print scenarionames
    print ' '
    print 'The comparison number for this comparison was: '
    print comparison_number
    print ' '
    print 'The resolutions of each of the scenarios were: '
    print scenarioresolution
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

    fg_transects_todo_file=[]
    num_scenarios = len(scenarionames)
    for j in range(num_scenarios):
        if (set(list(fg_types_todo[0])) == set(list(fg_types_todo[j]))): #scenarios same types
            for itype in fg_types_todo[0]:
                if (itype == 1):
                    fg_transects_todo_file.append(scenariodir[j] + '/fg_transects_todo.txt')
        else:
            print 'can not compare these scenarios. Their todo files have different types'
            print 'scenario number j different type from scenario number 0, j= ',j
            print 'fg_types_todo[0] was: ',fg_types_todo[0]
            print 'fg_types_todo[j] was: ',fg_types_todo[j]
            sys.exit(1)

    print 'the fg_transects_todo_file list was: '
    print fg_transects_todo_file
    print ' '

    ######  Now loop over all the work that has been done for these scenarios, and ###
    ######  check which transects to compare                                       ###
    transect_found=0
    for itype in fg_types_todo[0]:    #set(list(fg_types_todo[j]))=set(list(fg_types_todo[0]))##
        if (itype == 1):              #Do all the transect type of fixed grids
            transect_found=1
            fg_transects_todo=[]      #A list of arrays
            for j in range(num_scenarios):
                fg_transects_todo.append(loadtxt(fg_transects_todo_file[j],dtype=integer))
                if (size(fg_transects_todo[j]) == 1):
                    fg_transects_todo[j] = array([fg_transects_todo[j]])
                if (set(list(fg_transects_todo[0])) != set(list(fg_transects_todo[j]))):
                    print 'can not compare transects of these scenarios. '
                    print 'because they do not have the same transect numbers'
                    print 'fg_transects_todo[0] was: ',fg_transects_todo[0]
                    print 'j and fg_transects_todo[j] were: ',j,fg_transects_todo[j]
                    sys.exit(1)

            #All transects have the same numbers, so continue
            print 'the fg_transects_todo[0] was: '
            print fg_transects_todo[0]
            print ' '

            for iwhich in fg_transects_todo[0]:       #all fg_transects_todo are equal as a set
                print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
                print '++++++++++++ PROCESSING Transect number (iwhich) +++++++ ',iwhich
                print ' '
                hazard_curves_file=[]; OutZetafile=[]; #transect iwhich, multiple scenarios
                fixed_grid_file=[]

                for j in range(num_scenarios):
                    outdir=scenariodir[j]+'/'+scenarionames[j]+'_transect_fg'+str(iwhich)
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

                #The plotter will automatically choose evenly distributed points for
                #the hazard curves.  If these are good enough, choose whichrows=None;
                #otherwise, make whichrows a list (say of gauge indices) and the plotter
                #will append the evenly distributed points to the whichrows list.
                #BUT NOTE THAT the gauge indices can be transect dependent, and if this
                #is a priority, we would read them in from right scenario directory,
                #and make_scenario.py would have to be modified to put them there.

                whichrows=None
                depth=True; speed=False; mflux=False;
                zeta = linspace(0,12,1201)
                whichzetas=[0,100,200,300,400,600,800,1000]  #zeta=0.0,1.0,2.0,3.0,4.0,6.0,8.0,10.0
                whichzetas=array(whichzetas)
                pbarlist=[.1,.2,.3,.4,.5,.6,.7,.8,.9,.95]   #10 zeta-contours or 10 zeta-transects

                comparisondir=projectdir + '/comparisons/comparison_' + str(comparison_number) +\
                              '/transect_fg' + str(iwhich)
                print 'The comparison outputs for this transect will be written at: '
                print comparisondir
                print ' '
                print 'CALLING hazard_curves_transect_compare for fg transect: ',iwhich
                print ' '

                if 1:
                    hazard_curves_transect_compare(hazard_curves_file,fixed_grid_file,\
                          OutZetafile,\
                          zeta,whichzetas,comparisondir,pbarlist,\
                          CCtopo.X,CCtopo.Y,CCtopo.Z,scenarionames,scenarioresolution,\
                          whichrows=whichrows,speed=speed,depth=depth,mflux=mflux)

    if (transect_found == 0):
        print 'No transects were found to compare'
