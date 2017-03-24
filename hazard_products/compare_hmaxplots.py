##### Comparing the hmax data in all the .npy files in the run_num_list of the compare_dir_list to ###
##### the corresponding fine grid hmax data in the .npy file in the to_dir list.  Note that all    ###
##### or part of the runs for a given Mw category can be compared with one execution of this pgm   ###
##### but if more Mw categories are to be done, the USER INPUT below has to be changed and the     ###
##### pgm executed again.                                                                          ### 
#####
##### One .png plot (3 rows and 3 columns of sub plots) is stored in the run directory of the last ###
##### entry in the compare_dir_list (in this case, in the run directories for pseudoCM_40_9.2_runs)###
##### Each row gives the fine grid plot, the compare plot, and the absolute difference of the two. ###
##### Also left in this same directory is an output comments file called output.compare_hmaxplots. ### 
#####
##### To run:  python compare_hmaxplots.py     (run from $FEMA/ptha_tools directory)
#####                                          (when compare_hmaxplots.py is stored there)
from numpy import *
import pylab
from pylab import *
import os,sys
projectdir = os.environ['FEMA']

####################  USER INPUT                                            ############
####################                                                        ############
#################### SPECIFY WHICH run numbers to compare their hmax values ############
#################### Comparing the coarse, modified coarse, and pseudo hmax to fine ####
#################### only on the 2D grid, not the transects for now.                   #
########################################################################################
run_num_list = []
for j in range(100):
    run_num_list.append(j)
#run_num_list = [5]

#### Set which Mw is being created in this job run
Mw=str(9.2)

#### Set the compare_dir_list and compare_title_list
compare_title_list=[]; compare_dir_list=[];
compare_title_list.append('coarse')
compare_dir_list.append(projectdir + '/redclaw_1106/all_runs_npy_files/coarse_' + Mw + '_runs')
compare_title_list.append('coarse_mod')
compare_dir_list.append(projectdir + '/redclaw_1106/all_runs_npy_files/coarse_mod_' + Mw +'_runs')
compare_title_list.append('pseudoCM_40')
compare_dir_list.append(projectdir + '/redclaw_1106/all_runs_npy_files/pseudoCM_40_' + Mw +'_runs')
####
####  Always compare to fine, and use the fixed_grid_dir below
Ctitle = 'fine'
to_dir = projectdir + '/redclaw_1106/all_runs_npy_files/fine_' + Mw + '_runs'
fixed_grid_dir = projectdir + '/redclaw_1106/all_runs_npy_files'

####  Set the output comments file
comments = 'output.compare_hmaxplots'
####
#### END OF USER INPUT ###############################

### Dimensions for any 2Dgrids and aspect ratio for Crescent City
Npts1=88; Npts2=78;
aspectNO=41.75
###

fixed_big_grid_file = projectdir + '/DataFiles/' + 'CCtopo_1_3sec.tt3'
from clawpack.geoclaw.topotools import Topography
CCtopo = Topography()
CCtopo.read(fixed_big_grid_file,3)

fixed_grid_file = fixed_grid_dir + '/fine_xyB0_fg3.npy'
d=load(fixed_grid_file)
x=d[:,0]; y=d[:,1];
xsquare=x.reshape(Npts1,Npts2,order='F')
ysquare=y.reshape(Npts1,Npts2,order='F')

#fixed_grid_fileC = fixed_grid_dir + '/coarse_xyB0_fg3.npy' #do not need, not using eta
#B0_fine=d[:,2];                                            #do not need, since not finding eta
#dc=load(fixed_grid_fileC)                                  #x and y same from fine and coarse
#B0_coarse=dc[:,2];                                         #do not need, since not finding eta

region_extent=[]
region_extent.append(x[0]); region_extent.append(x[-1]);
region_extent.append(y[0]); region_extent.append(y[-1]);
    
figure(12,(12,8))
print ' '
print 'The output comments for this job run are in the run directories in the directory: '
print compare_dir_list[-1]
print ' '
for run_num in run_num_list:
    commentsfile = compare_dir_list[-1] + '/run_' + str(run_num) + '/' + comments
    commentsout  = open(commentsfile,'w')
    commentsout.write(' \n')
    commentsout.write('   HMAX COMPARISONS OF RUN_%s TO %s RUN_%s  \n' %(run_num,Ctitle,run_num))
    commentsout.write(' \n')
    commentsout.write(' \n')
    tostr = to_dir + '/run_' + str(run_num)
    commentsout.write('To directory:  \n')
    commentsout.write('    %s \n' %tostr)
    commentsout.write(' \n')

    commentsout.write('Comparing hmax in directories below to fine above:  \n')
    for k in [0,1,2]:
        comparestr = compare_dir_list[k] + '/run_' + str(run_num)
        commentsout.write('    %s \n' %comparestr)
        commentsout.write(' \n')
    to_file = to_dir + '/run_' + str(run_num) + '/hmax_fg3.npy'
    out_file = compare_dir_list[-1] +'/run_'+str(run_num)+'/hmax_diff_to_' + Ctitle +'.png' 

    clf()
    to_data=load(to_file)
    cmaxto=max(to_data)
    cmeanto=mean(to_data)
    to_data_square=to_data.reshape(Npts1,Npts2,order='F')

    plotting_no = 0
    for k in [0,1,2]:
        ## Load the compare file
        compare_file = compare_dir_list[k] + '/run_' + str(run_num) + '/hmax_fg3.npy'
        compare_title=compare_title_list[k]
        compare=load(compare_file)
        diff=abs(to_data-compare)
        cmaxdiff=max(diff)
        cmeandiff=mean(diff)
        cmaxcompare=max(compare)
        cmeancompare=mean(compare)
        diff_square=diff.reshape(Npts1,Npts2,order='F')
        compare_square=compare.reshape(Npts1,Npts2,order='F')

        commentsout.write('  abs(%s-%s): max and mean were %s %s  \n' %(Ctitle,compare_title,cmaxdiff,cmeandiff))
        commentsout.write('  %s: hmax_max, hmax_mean were %s %s  \n' %(Ctitle,cmaxto,cmeanto))
        commentsout.write('  %s: hmax_max, hmax_mean were %s %s  \n' %(compare_title,cmaxcompare,cmeancompare))
        commentsout.write(' \n')
        if (cmaxdiff <= 16.0):
            clinesdiff = [1e-3] + [0.5,1.0] + list(linspace(2.0,16.0,8))
        else:
            clinesdiff = [1e-3] + [0.5,1.0] + list(linspace(2.0,14.0,7)) + [cmaxdiff]
        if (cmaxto <= 16.0):
            clinesto = [1e-3] + [0.5,1.0] + list(linspace(2.0,16.0,8))
        else:
            clinesto = [1e-3] + [0.5,1.0] + list(linspace(2.0,14.0,7)) + [cmaxto]
        if (cmaxcompare <= 7.0):
            clinescompare = [1e-3] + [0.5,1.0] + list(linspace(2.0,16.0,8))
        else:
            clinescompare = [1e-3] + [0.5,1.0] + list(linspace(2.0,14.0,7)) + [cmaxcompare]

        nlines = len(clinesto)
        n1 = int(floor((nlines-1)/2.))
        n2 = nlines - 1 - n1
        Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
        Red = hstack([linspace(0,0.8,n1), ones(n2)])
        Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
        colors2 = zip(Red,Green,Blue)

        #### Do the to_data subplot
        plotting_no += 1
        subplot(3,3,plotting_no)
        CS100=contourf(xsquare,ysquare,to_data_square,\
                       clinesto,colors=colors2)
        cbar_shrink=1.0
        cbar = colorbar(shrink=cbar_shrink)
        cbar.set_ticks(clinesto)
        hold(True)

        ##Contours of topo added:
        if 1:
            #contour(CCtopo.X,CCtopo.Y,CCtopo.Z, arange(0,21,2), colors='k')
            contour(CCtopo.X,CCtopo.Y,CCtopo.Z, arange(0,8,4), colors='k')
            #xlim(-124.215,-124.18)
            xlim(region_extent[0],region_extent[1])
            if (plotting_no > 6):
                xticks(rotation=20)
            else:
                xticks([])
            #ylim(41.735,41.768)
            ylim(region_extent[2],region_extent[3])
            ticklabel_format(format='plain',useOffset=False)
        a = gca()
        a.set_aspect(1./cos(aspectNO*pi/180.))
        title('hmax(%s), run=%s' %(Ctitle,run_num))

        ### Do the compare subplot
        plotting_no += 1
        subplot(3,3,plotting_no)
        CS100=contourf(xsquare,ysquare,compare_square,\
                       clinescompare,colors=colors2)
        cbar_shrink=1.0
        cbar = colorbar(shrink=cbar_shrink)
        cbar.set_ticks(clinescompare)
        hold(True)

        ##Contours of topo added:
        if 1:
            #contour(CCtopo.X,CCtopo.Y,CCtopo.Z, arange(0,21,2), colors='k')
            contour(CCtopo.X,CCtopo.Y,CCtopo.Z, arange(0,8,4), colors='k')
            #xlim(-124.215,-124.18)
            xlim(region_extent[0],region_extent[1])
            if (plotting_no > 6):
                xticks(rotation=20)
            else:
                xticks([])
            #ylim(41.735,41.768)
            ylim(region_extent[2],region_extent[3])
            ticklabel_format(format='plain',useOffset=False)
        a = gca()
        a.set_aspect(1./cos(aspectNO*pi/180.))
        title('hmax(%s), run=%s' %(compare_title,run_num))

        ### Do the diff subplot
        plotting_no += 1
        subplot(3,3,plotting_no)
        CS100=contourf(xsquare,ysquare,diff_square,\
                       clinesdiff,colors=colors2)
        cbar_shrink=1.0
        cbar = colorbar(shrink=cbar_shrink)
        cbar.set_ticks(clinesdiff)
        hold(True)

        ##Contours of topo added:
        if 1:
            #contour(CCtopo.X,CCtopo.Y,CCtopo.Z, arange(0,21,2), colors='k')
            contour(CCtopo.X,CCtopo.Y,CCtopo.Z, arange(0,8,4), colors='k')
            #xlim(-124.215,-124.18)
            xlim(region_extent[0],region_extent[1])
            if (plotting_no > 6):
                xticks(rotation=20)
            else:
                xticks([])
            #ylim(41.735,41.768)
            ylim(region_extent[2],region_extent[3])
            ticklabel_format(format='plain',useOffset=False)
        a = gca()
        a.set_aspect(1./cos(aspectNO*pi/180.))
        #title('abs(hmax(%s) - hmax(coarse_mod)), run=%s, fg=%s' %(Ctitle,run_num,fg))
        title('abs(h diffs), run=%s' %run_num)

    savefig(out_file)

