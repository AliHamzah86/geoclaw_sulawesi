#################################################
#### This program creates the event probs
#### for the clustered data stored in the cluster_file
#### with associated master probs stored in the
#### master_prob_wts_file. The output is a text file
#### of the event probabilities that are needed by
#### the make_scenario.py program.
#### 
#### Input:  cluster_file
####         master_prob_wts_file
####         num_clusters  (number of clusters)
####         mag_factors   (dictionary of factors used for each Mw category, they sum to 1)
####
#### Output (example): output.cluster_event_probs_make_uniform_etaClustering1106 
####
#### These event probs are used in the make_scenario.py program.
#### At the moment, they are manually put into this program
#### which will then store them in the proper scenario directory's
#### scenario_prb_wgts.txt file in the proper format for further use.
####
#### To run:
####  python cluster_event_probs_make.py > output.cluster_event_probs_make_uniform_etaClustering1106
####

from numpy import *
import os,sys
projectdir = os.environ['FEMA']

###############    USER INPUT ##################
mag_factors={}
#num_clusters = 20
num_clusters = 40

if 1:  #new, 400
    mag_factors[9.2]=.1;  mag_factors[9.0]=.3;
    mag_factors[8.8]=.3; mag_factors[8.6]=.3;
    print 'Calculating scenario probs for clustering with the 400 files'
    print ' '

if 0:  #old, 100, all were 9.0 magnitudes
    mag_factors[9.0]=1.0 
    print 'Calculating scenario probs for clustering with the 100 files'
    print ' '

if 1:
    #cluster_file = projectdir + '/scenarios_new/scenario_dzCC_dtopo_etamax_20cls/clusters.txt'
    #cluster_file = projectdir + '/scenarios_1106/scenario_dzCC_dtopo_etamax_20cls/clusters.txt'
    #cluster_file = projectdir + '/scenarios_1106/scenario_dtopo_40cls/clusters.txt'
    #cluster_file = projectdir + '/scenarios_1106/scenario_all_pseudoCM_dzCC_PE/clusters.txt'
    cluster_file = projectdir + '/scenarios_1106/scenario_all_pseudoCM_dzCC_PE_40/clusters.txt'
    #cluster_file = projectdir + '/scenarios_1106/scenario_all_pseudoCM/clusters.txt'
    #cluster_file = projectdir + '/scenarios_1106/scenario_all_pseudoCM_40/clusters.txt'

if 0:
    #cluster_file = projectdir + '/scenarios_new/scenario_coarse_eta_20cls/clusters.txt'
    #cluster_file = projectdir + '/scenarios_1106/scenario_coarse_eta_20cls/clusters.txt'
    cluster_file = projectdir + '/scenarios_1106/scenario_coarse_eta_40cls/clusters.txt'

if 1: #not uniform - what Don Sub made
    print 'The master prob wts file was Don Sub made, not uniform'
    master_prob_wts_file = projectdir + '/scenario_prb_wgts.txt'

if 0: #uniform within a grouping of 100 of same magnitude
    print 'The master prob wts file was uniform'
    master_prob_wts_file = projectdir + '/scenario_prb_wgts_uniform.txt'

####  END OF USER INPUT #####################

print ' '
print 'The cluster_file used was: '
print cluster_file
print ' '
print 'The master_prob_wts file was: '
print master_prob_wts_file
print ' '

## Load the clusters_data
clusters_data = loadtxt(cluster_file,dtype=integer)

## Load the correct master_prob_wts
master_prob_wts = loadtxt(master_prob_wts_file)

j=0
final_run_numbers=[]; final_mags=[]; final_probs=[];
for ncls in range(num_clusters):
    cluster_no = clusters_data[j,0]
    number_in_cluster = clusters_data[j,1]
    print ' '
    print 'PROCESSING cluster_no, number_in_cluster: ',cluster_no, number_in_cluster
    j +=1
    centroid_run_number = clusters_data[j,0]
    centroid_mag = clusters_data[j,1]/10.0
    final_run_numbers.append(centroid_run_number)
    final_mags.append(clusters_data[j,1])

    j +=1
    probsum=0.0
    for imember in range(number_in_cluster):
        member_run_number = clusters_data[j,0];
        member_mag = clusters_data[j,1]/10.0
        print '   Processing member, mag: ',member_run_number,member_mag
        probsum += mag_factors[member_mag]*master_prob_wts[member_run_number]
        j +=1
    print '       Cluster Centroid Run Number: ',centroid_run_number
    print '       Cluster Centroid Magnitude: ',centroid_mag
    print '       Cluster probability: ',probsum
    final_probs.append(probsum)

print ' '
print 'The final_run_numbers were: '
print final_run_numbers
print ' '
print 'The magnitudes*10 for the cluster centroids were: '
print final_mags
print ' '
print 'The final scenario probability weights were: '
print final_probs
print ' '
print 'The number of lines read in the cluster file was: ',j
