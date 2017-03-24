##### Compares the .npy files in the run_num_list in the directory compare_dir to the           ####
##### corresponding .npy file in the directory to_dir. The to_dir  can be for either a fine     ####
##### grid run or a coarse grid run.  If fine, Ctitle='fine'; otherwise, Ctitle='coarse'.       ####
##### The title for the plots for the compare file is set in compare_title.                     ####
##### The program is presently set up to compare fine_8.6_runs to the coarse_8.6_runs           ####
##### just for run_5.                                                                           ####
#####
##### A comments file called output.compare_to_coarse is left in the compare_dir/run_5.         ####
##### If Ctitle were 'fine' (say we were comparing coarse_mod to fine), then the comments       ####
##### file would be called output.compare_to_fine.  Can only compare "to" coarse or fine.       ####
##### Three .png files are also stored in to_dir/run_5 for this example.  They are for          ####
##### the comparisons on Transect 1, Transect 2, and the 2D grid and are called                 ####
##### hmax_diff_to_coarse_fg1.png, hmax_diff_to_coarse_fg2.png, and hmax_diff_to_coarse.fg3.png ####
#####
##### To run:  python compare_to.py     (run from $FEMA/ptha_tools directory)
#####                                   (when compare_to.py is stored there)

from numpy import *
import pylab
from pylab import *
import os,sys
projectdir = os.environ['FEMA']

#####            USER INPUT                                                 ###########
#####                                                                       ###########
#################### SPECIFY WHICH run numbers to compare their hmax values ###########
####################  comparing the to_dir run_num to the compare_dir run number ####
####################                                                                  #
#######################################################################################
run_num_list = []
#for j in range(100):
#    run_num_list.append(j)
run_num_list = [5]

if 0:
    #Those without the 40 extension were made with 20 clusters
    #compare_title = 'pseudo'
    #compare_title = 'pseudoCM'
    compare_title = 'pseudoCM_dzCC_PE'
    #compare_title = 'pseudo_40'
    #compare_dir = projectdir + '/all_runs_npy_files_new/pseudo_9.2_runs'
    #compare_dir = projectdir + '/redclaw_1106/all_runs_npy_files/pseudo_9.2_runs'
    compare_dir = projectdir + '/redclaw_1106/all_runs_npy_files/pseudoCM_dzCC_PE_9.2_runs'
    #compare_dir = projectdir + '/redclaw_1106/all_runs_npy_files/pseudo_40_9.2_runs'
if 0:
    compare_title = 'pseudo_dtopo'
    #compare_dir = projectdir + '/all_runs_npy_files_new/pseudo_dtopo_9.2_runs'
    compare_dir = projectdir + '/redclaw_1106/all_runs_npy_files/pseudo_dtopo_9.2_runs'
if 0:
    compare_title = 'coarse_mod'
    #compare_dir = projectdir + '/all_runs_npy_files_new/coarse_mod_9.2_runs'
    compare_dir = projectdir + '/redclaw_1106/all_runs_npy_files/coarse_mod_9.2_runs'

if 1:
    compare_title = 'fine'
    compare_dir = projectdir + '/redclaw_1106/all_runs_npy_files/fine_8.6_runs'

#fixed_grid_dir = projectdir + '/all_runs_npy_files_new'
fixed_grid_dir = projectdir + '/redclaw_1106/all_runs_npy_files'

## Choose which we are comparing to, either coarse or fine runs
if 1:  #which we are comparing to, either coarse or fine
    Ctitle = 'coarse'
    #to_dir = projectdir + '/all_runs_npy_files_new/coarse_9.2_runs'
    to_dir = projectdir + '/redclaw_1106/all_runs_npy_files/coarse_8.6_runs'
    comments = 'output.compare_to_coarse'
if 0:
    Ctitle = 'fine'
    #to_dir = projectdir + '/all_runs_npy_files_new/fine_9.2_runs'
    to_dir = projectdir + '/redclaw_1106/all_runs_npy_files/fine_9.2_runs'
    comments = 'output.compare_to_fine'
#####     END OF USER INPUT ##########################################################

### Dimensions for any 2Dgrids and aspect ratio for Crescent City
Npts1=88; Npts2=78;
aspectNO=41.75
###

fixed_big_grid_file = projectdir + '/DataFiles/' + 'CCtopo_1_3sec.tt3'
from clawpack.geoclaw.topotools import Topography
CCtopo = Topography()
CCtopo.read(fixed_big_grid_file,3)

to_file={}; compare_file={}; out_file={};
x={}; y={}; B0_fine={}; B0_coarse={};   #may need B0_fine, B0_coarse for eta

for key in [1,2,3]:
    fixed_grid_file = fixed_grid_dir + '/fine_xyB0_fg' + str(key) + '.npy'
    fixed_grid_fileC = fixed_grid_dir + '/coarse_xyB0_fg' + str(key) + '.npy'
    d=load(fixed_grid_file)
    x[key]=d[:,0]; y[key]=d[:,1]; B0_fine[key]=d[:,2];
    dc=load(fixed_grid_fileC)          #x and y same from fine and coarse files
    B0_coarse[key]=dc[:,2];

region_extent=[]
region_extent.append(x[3][0]); region_extent.append(x[3][-1]);
region_extent.append(y[3][0]); region_extent.append(y[3][-1]);
    
figure(12,(12,8))
for run_num in run_num_list:
    commentsfile = compare_dir + '/run_' + str(run_num) + '/' + comments
    commentsout  = open(commentsfile,'w')
    commentsout.write(' \n')
    commentsout.write('   COMPARISON OF RUN_%s TO %s RUN_%s  \n' %(run_num,Ctitle,run_num))
    commentsout.write(' \n')
    commentsout.write(' \n')
    comparestr = compare_dir + '/run_' + str(run_num)
    tostr = to_dir + '/run_' + str(run_num)
    commentsout.write('To directory:  \n')
    commentsout.write('    %s \n' %tostr)
    commentsout.write(' \n')
    commentsout.write('Compare and diff plots directory:  \n')
    commentsout.write('    %s \n' %comparestr)
    commentsout.write(' \n')
    for key in [1,2,3]:
        to_file[key] = to_dir + '/run_' + str(run_num) + '/hmax_fg' + str(key) + '.npy'
        compare_file[key] = compare_dir + '/run_' + str(run_num) + '/hmax_fg' + str(key) + '.npy'
        out_file[key] = compare_dir +'/run_'+str(run_num)+'/hmax_diff_to_' + Ctitle +'_fg'+str(key)+'.png' 

    for fg in [1,2]:
        clf()
        to_data=load(to_file[fg])
        if (Ctitle == 'coarse'):
            todata_eta = to_data + B0_coarse[fg]
        else:
            todata_eta = to_data + B0_fine[fg]
        compare=load(compare_file[fg])
        compare_eta = compare + B0_fine[fg]
        diff=abs(to_data-compare)
        diffmax=max(diff); diffmean=mean(diff);
        diff_eta=abs(todata_eta-compare_eta)
        diffmax_eta=max(diff_eta); diffmean_eta=mean(diff_eta);
        commentsout.write('Processing fixed grid no = %s  \n' %fg)
        commentsout.write('  diffmax and diffmean were %s %s  \n' %(diffmax,diffmean))
        commentsout.write('  diffmax_eta and diffmean_eta were %s %s  \n ' %(diffmax_eta,diffmean_eta))
        commentsout.write(' \n')
        yplot=y[fg]
        plot(yplot,diff,'r')
        plot(yplot,diff_eta,'b')
        #xticks([yplot[0],yplot[-1]],[str(yplot[0]),str(yplot[-1])],rotation=20)
        xticks(rotation=20)
        xlabel('Latitude on the Transect')
        ylabel('abs(differences) in meters: Red(h), Blue(eta)')
        title('abs differences (%s - %s): Red(h), Blue(eta), run=%s, fg=%s'\
                                %(compare_title,Ctitle,run_num,fg))
        savefig(out_file[fg])
    
    clf()
    fg=3
    to_data=load(to_file[fg])
    compare=load(compare_file[fg])

    ######  Compute eta ##########
    compare_eta = compare + B0_fine[fg]
    if (Ctitle == 'coarse'):
        todata_eta = to_data + B0_coarse[fg]
    else:
        todata_eta = to_data + B0_fine[fg]
    ##############################

    diff=abs(to_data-compare)
    diff_eta=abs(todata_eta-compare_eta)

    cmaxdiff=max(diff)
    cmeandiff=mean(diff)
    cmaxdiff_eta=max(diff_eta)
    cmeandiff_eta=mean(diff_eta)
    cmaxto=max(to_data)
    cmeanto=mean(to_data)
    cmaxcompare=max(compare)
    cmeancompare=mean(compare)

    diff=diff.reshape(Npts1,Npts2,order='F')
    diff_eta=diff_eta.reshape(Npts1,Npts2,order='F')
    to_data=to_data.reshape(Npts1,Npts2,order='F')
    compare=compare.reshape(Npts1,Npts2,order='F')
    xsquare=x[fg].reshape(Npts1,Npts2,order='F')
    ysquare=y[fg].reshape(Npts1,Npts2,order='F')

    commentsout.write('Processing 2Dgrid no = %s \n' %fg)
    commentsout.write('  abs(%s-compare): max and mean were %s %s  \n' %(Ctitle,cmaxdiff,cmeandiff))
    commentsout.write('  abs(%s_eta-compare_eta): max and mean were %s %s  \n'\
                         %(Ctitle,cmaxdiff_eta,cmeandiff_eta))
    commentsout.write('  %s: hmax_max, hmax_mean were %s %s  \n' %(Ctitle,cmaxto,cmeanto))
    commentsout.write('  compare: hmax_max, hmax_mean were %s %s  \n' %(cmaxcompare,cmeancompare))
    commentsout.write(' \n')
    if (cmaxdiff <= 16.0):
        clinesdiff = [1e-3] + [0.5,1.0] + list(linspace(2.0,16.0,8))
    else:
        clinesdiff = [1e-3] + [0.5,1.0] + list(linspace(2.0,14.0,7)) + [cmaxdiff]
    if (cmaxdiff_eta <= 16.0):
        clinesdiff_eta = [1e-3] + [0.5,1.0] + list(linspace(2.0,16.0,8))
    else:
        clinesdiff_eta = [1e-3] + [0.5,1.0] + list(linspace(2.0,14.0,7)) + [cmaxdiff_eta]
    if (cmaxto <= 16.0):
        clinesto = [1e-3] + [0.5,1.0] + list(linspace(2.0,16.0,8))
    else:
        clinesto = [1e-3] + [0.5,1.0] + list(linspace(2.0,14.0,7)) + [cmaxto]
    if (cmaxcompare <= 7.0):
        clinescompare = [1e-3] + [0.5,1.0] + list(linspace(2.0,16.0,8))
    else:
        clinescompare = [1e-3] + [0.5,1.0] + list(linspace(2.0,14.0,7)) + [cmaxcompare]

    ### Do the diff subplot
    nlines = len(clinesdiff)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
    Red = hstack([linspace(0,0.8,n1), ones(n2)])
    Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
    colors2 = zip(Red,Green,Blue)

    subplot(2,2,1)
    CS100=contourf(xsquare,ysquare,diff,\
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
        xticks([])
        #ylim(41.735,41.768)
        ylim(region_extent[2],region_extent[3])
        ticklabel_format(format='plain',useOffset=False)
    a = gca()
    a.set_aspect(1./cos(aspectNO*pi/180.))
    #title('abs(hmax(%s) - hmax(coarse_mod)), run=%s, fg=%s' %(Ctitle,run_num,fg))
    title('abs(h diffs), run=%s, fg=%s' %(run_num,fg))

    ### Do the diff_eta subplot
    nlines = len(clinesdiff_eta)
    n1 = int(floor((nlines-1)/2.))
    n2 = nlines - 1 - n1
    Green = hstack([linspace(1,1,n1),linspace(1,0,n2)])
    Red = hstack([linspace(0,0.8,n1), ones(n2)])
    Blue = hstack([linspace(1,0.2,n1), zeros(n2)])
    colors2 = zip(Red,Green,Blue)

    subplot(2,2,2)
    CS100=contourf(xsquare,ysquare,diff_eta,\
                   clinesdiff_eta,colors=colors2)
    cbar_shrink=1.0
    cbar = colorbar(shrink=cbar_shrink)
    cbar.set_ticks(clinesdiff_eta)
    hold(True)

    ##Contours of topo added:
    if 1:
        #contour(CCtopo.X,CCtopo.Y,CCtopo.Z, arange(0,21,2), colors='k')
        contour(CCtopo.X,CCtopo.Y,CCtopo.Z, arange(0,8,4), colors='k')
        #xlim(-124.215,-124.18)
        xlim(region_extent[0],region_extent[1])
        xticks([])
        #ylim(41.735,41.768)
        ylim(region_extent[2],region_extent[3])
        ticklabel_format(format='plain',useOffset=False)
    a = gca()
    a.set_aspect(1./cos(aspectNO*pi/180.))
    title('abs(eta diffs), run=%s, fg=%s' %(run_num,fg))


    #### Do the to_data subplot
    subplot(2,2,3)
    CS100=contourf(xsquare,ysquare,to_data,\
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
        xticks(rotation=20)
        #ylim(41.735,41.768)
        ylim(region_extent[2],region_extent[3])
        ticklabel_format(format='plain',useOffset=False)
    a = gca()
    a.set_aspect(1./cos(aspectNO*pi/180.))
    title('hmax(%s), run=%s, fg=%s' %(Ctitle,run_num,fg))

    ### Do the compare subplot
    subplot(2,2,4)
    CS100=contourf(xsquare,ysquare,compare,\
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
        xticks(rotation=20)
        #ylim(41.735,41.768)
        ylim(region_extent[2],region_extent[3])
        ticklabel_format(format='plain',useOffset=False)
    a = gca()
    a.set_aspect(1./cos(aspectNO*pi/180.))
    title('hmax(%s), run=%s, fg=%s' %(compare_title,run_num,fg))
    savefig(out_file[fg])

