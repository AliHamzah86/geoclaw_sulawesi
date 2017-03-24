#### This program turns the 400 coarse_mod .npy files 
#### into the 400 pseudo fine .npy files. For the 20 (or 40) centroids
#### the pseudo files will be the fine .npy files.
#### In the USER INPUT, the cluster_file used to do this is specified,
#### as well as the pseudoname indicating which directory to put the
#### output pseudo .npy file. The number of clusters used in the cluster
#### file is also specified. Note that all the Mw categories are included
#### in the cluster_file, so all 400 are done.
#### 
#### A comments file called output.pseudoCM_dzCC_PE_40_create gives details
#### of the job run and some comparison data.  (If pseudoname were pseudoCM_40,
#### we would store the comments on output.pseudoCM_40_create for example.)
####
#### To run:   Go to directory $FEMA/redclaw_1106/all_runs_npy_files and do:
####           python $FEMA/ptha_tools/pseudo_create.py > output.pseudoCM_dzCC_PE_40_create
####
####           when all python codes have been stored in $FEMA/ptha_tools.
################################################################################
from numpy import *
import os,sys

projectdir = os.environ['FEMA']
#runsdir = 'all_runs_npy_files_new'
runsdir = 'redclaw_1106/all_runs_npy_files'  #always choose for the FEMA 1106 results

###############  USER INPUT   #########################################################
if 0: #uses the etamean-etamax clustering of the coarse_mod (CM) .npy files rather than the coarse files
    #cluster_file = projectdir + '/scenarios_1106/scenario_all_pseudoCM/clusters.txt'
    #pseudoname = 'pseudoCM_'
    cluster_file = projectdir + '/scenarios_1106/scenario_all_pseudoCM_40/clusters.txt'
    pseudoname = 'pseudoCM_40_'
if 1: #uses the new dzCC-PE clustering of the coarse_mod (CM) .npy files
    #cluster_file = projectdir + '/scenarios_1106/scenario_all_pseudoCM_dzCC_PE/clusters.txt'
    #pseudoname = 'pseudoCM_dzCC_PE_'
    cluster_file = projectdir + '/scenarios_1106/scenario_all_pseudoCM_dzCC_PE_40/clusters.txt'
    pseudoname = 'pseudoCM_dzCC_PE_40_'
if 0: #uses the etamean-etamax clustering
    #cluster_file = projectdir + '/scenarios_new/scenario_coarse_eta_20cls/clusters.txt'
    #cluster_file = projectdir + '/scenarios_1106/scenario_coarse_eta_20cls/clusters.txt'
    #pseudoname = 'pseudo_'
    cluster_file = projectdir + '/scenarios_1106/scenario_coarse_eta_40cls/clusters.txt'
    pseudoname = 'pseudo_40_'

if 0: #uses the dtopo clustering
    #cluster_file = projectdir + '/scenarios_new/scenario_dzCC_dtopo_etamax_20cls/clusters.txt'
    #cluster_file = projectdir + '/scenarios_1106/scenario_dzCC_dtopo_etamax_20cls/clusters.txt'
    #pseudoname = 'pseudo_dtopo_'
    cluster_file = projectdir + '/scenarios_1106/scenario_dtopo_40cls/clusters.txt'
    pseudoname = 'pseudo_dtopo_40_'

#### Set the number of clusters
#num_clusters = 20     #number of clusters
num_clusters = 40      #number of clusters
######                 END OF USER INPUT        ##############

print 'The cluster_file used was: '
print cluster_file
print ' '
print 'The directory for the fine, coarse_mod, pseudo .npy files is: '
dirstr=projectdir + '/' + runsdir
print dirstr
print ' '

## Load the clusters_data
clusters_data = loadtxt(cluster_file,dtype=integer)

j=0
for ncls in range(num_clusters):
    cluster_no = clusters_data[j,0]
    number_in_cluster = clusters_data[j,1]
    print 'PROCESSING cluster_no, number_in_cluster: ',cluster_no, number_in_cluster
    j +=1
    centroid_run_number = clusters_data[j,0]
    centroid_mag = clusters_data[j,1]/10.0
    print '   centroid_no, centroid_mag were: ',centroid_run_number, centroid_mag

    centroid_fine_dir = projectdir + '/' + runsdir + '/fine_' + str(centroid_mag) +\
                               '_runs/run_' + str(centroid_run_number)
    centroid_coarse_mod_dir = projectdir + '/' + runsdir + '/coarse_mod_' + str(centroid_mag) +\
                               '_runs/run_' + str(centroid_run_number)
    if 0:
        print '   centroid fine_dir: '
        print '   ',centroid_fine_dir
        print '   centroid_coarse_mod_dir: '
        print '   ',centroid_coarse_mod_dir
    delta={};
    for fg in [1,2,3]:
        hmax_cen_fine=load(centroid_fine_dir + '/hmax_fg' + str(fg) + '.npy')
        hmax_cen_coarse_mod=load(centroid_coarse_mod_dir + '/hmax_fg' + str(fg) + '.npy')
        delta[fg]=hmax_cen_fine - hmax_cen_coarse_mod

    j +=1
    for imember in range(number_in_cluster):
        member_run_number = clusters_data[j,0];
        member_mag = clusters_data[j,1]/10.0
        print '   Processing member, mag: ',member_run_number,member_mag
        member_coarse_mod_dir = projectdir + '/' + runsdir + '/coarse_mod_' + str(member_mag) +\
                               '_runs/run_' + str(member_run_number)
        if 0:
            print '    member coarse mod dir: '
            print '    ',member_coarse_mod_dir

        ### for comparison only
        member_fine_dir = projectdir + '/' + runsdir + '/fine_' + str(member_mag) +\
                               '_runs/run_' + str(member_run_number)
        if 0:
            print '    member fine dir for comparison: '
            print '    ',member_fine_dir
            print ' '
        ###

        pseudo_filedir= projectdir + '/' + runsdir + '/' + pseudoname + str(member_mag) +\
                         '_runs/run_' + str(member_run_number)
        for fg in [1,2,3]:
            hmax_member_coarse_mod=load(member_coarse_mod_dir + '/hmax_fg' + str(fg) + '.npy')
            hmax_fine=load(member_fine_dir + '/hmax_fg' + str(fg) + '.npy')

            ###  create pseudo fine fine here, don't let hmax_pseudo drop below zero
            hmax_pseudo = hmax_member_coarse_mod + delta[fg]
            hmax_pseudo = where(hmax_pseudo < 0.0, 0.0, hmax_pseudo)
            pseudo_file= pseudo_filedir + '/hmax_fg' + str(fg) + '.npy'
            save(pseudo_file,hmax_pseudo)
            ###

            if 1:
                print '      fg, max(abs(diff)) between fine and pseudo fine was: ',\
                             fg, max(abs(hmax_pseudo-hmax_fine))
                print '      fg, max(abs(diff)) between fine and coarse_mod was: ',\
                             fg, max(abs(hmax_member_coarse_mod-hmax_fine))
        if 1:
            print ' '
            print '      Saved pseudo .npy files in : '
            print '        ',pseudo_filedir
            print ' '
        print ' '
        j +=1

print 'The number of lines read in the cluster file was: ',j
