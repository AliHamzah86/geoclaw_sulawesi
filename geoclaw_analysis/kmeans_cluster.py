"""
k-means clustering based on cluster_vars, e.g.
dtopo_proxies or coarse_mod_eta
"""

import matplotlib
matplotlib.use('Agg')  # for non-interactive

from pylab import *

try:
    from sklearn.cluster import KMeans
except:
    print "You need to install sklearn first"

all_runs_dir = '../geoclaw_output/all_runs_npy_files'

run_probs = loadtxt('../geoclaw_driver/scenario_prb_wgts.txt')  # 100 weights for each runno

# Relative weight for each Magnitude:
Mw_weights = {8.6:0.3, 8.8:0.3, 9.0:0.3, 9.2:0.1}


#cluster_vars = 'dtopo_proxies'
#cluster_vars = 'coarse_eta'
#cluster_vars = 'dzCC_PE'
#cluster_vars = 'coarse_mod_eta'

cluster_vars_list = ['dtopo_proxies', 'coarse_eta']

save_png = True
make_files = True


Mwlist = [8.6,8.8,9.0,9.2]
Mw_colors = {8.6:'m', 8.8:'g', 9.0:'b', 9.2:'r'}

for cluster_vars in cluster_vars_list:

    resolution = 'coarse'
    if cluster_vars == 'coarse_mod_eta':
        resolution = 'coarse_mod'
    fgno = 3

    events_description = []
    events_description2 = []
    events_prob = []
    events_Mw = []
    events_runno = []

    # from GeoClaw runs:
    events_etamean = []
    events_etamax = []
    events_hmax_sum = []

    # for dtopo proxies:
    events_dzCC = []
    events_dtopo_etamax = []
    events_PE = []
    events_dV = []

    for Mw in Mwlist:

        fname = all_runs_dir + '/%s_%s_etamax_fg%s.txt' % (resolution, Mw, fgno)
        geoclaw_output = loadtxt(fname, skiprows=1)

        fname = all_runs_dir + '/dtopo_proxies_%s.txt' % Mw
        dtopo_proxies = loadtxt(fname, skiprows=1)

        for runno in range(100):
            events_runno.append(runno)
            events_Mw.append(Mw)

            # probability for this event:
            events_prob.append(run_probs[runno] * Mw_weights[Mw])

            events_description.append("Mw=%3.1f, run=%3i" % (Mw,runno))
            events_description2.append("%3i   %3i" % (runno, (Mw*10)))
            events_etamean.append(geoclaw_output[runno,2])
            events_etamax.append(geoclaw_output[runno,3])
            events_hmax_sum.append(geoclaw_output[runno,4])

            events_dzCC.append(dtopo_proxies[runno,1])
            events_dtopo_etamax.append(dtopo_proxies[runno,2])
            events_PE.append(dtopo_proxies[runno,3])
            events_dV.append(dtopo_proxies[runno,4])
        
    events_dtopo_etamax = array(events_dtopo_etamax)

    # k-means clustering based on coarse grid etamax and etamean:


    n_clusters = 20


    if cluster_vars=='coarse_eta':
        cluster_data = vstack((events_etamean,events_etamax)).T
        scenario_name = 'coarse_eta_%sclusters' % n_clusters
        v1label = 'etamean'
        v2label = 'etamax'
        plot_title = 'k-means clustering based on coarse grid proxies'
    elif cluster_vars=='dtopo_proxies':
        cluster_data = vstack((events_dzCC,events_dtopo_etamax/3.)).T
        scenario_name = 'dzCC_dtopo_etamax_%sclusters' % n_clusters
        v1label = 'dzCC'
        v2label = 'dtopo_etamax/3'
        plot_title = 'k-means clustering based on dtopo proxies'
    elif cluster_vars=='dzCC_PE':
        cluster_data = vstack((events_dzCC,log(array(events_PE)))).T
        scenario_name = 'dzCC_PE_%sclusters' % n_clusters
        v1label = 'dzCC'
        v2label = 'log PE'
        plot_title = 'k-means clustering based on dzCC and PE'
    elif cluster_vars=='coarse_mod_eta':
        cluster_data = vstack((events_etamean,events_etamax)).T
        scenario_name = 'coarse_mod_eta_%sclusters' % n_clusters
        v1label = 'etamean'
        v2label = 'etamax'
        plot_title = 'k-means clustering based on coarse mod proxies'
    else:
        raise InputError('*** Unrecognized cluster_vars = %s' % cluster_vars)

    #plot_title = 'k-means clustering based on %s and %s' % (v1label,v2label)


    # compute k-mean clusters:
    random_state = 12345  # fix seed for reproducibility
    kmeans = KMeans(n_clusters=n_clusters, max_iter=10000,\
            n_init=100,random_state=random_state).fit(cluster_data)


    centers = kmeans.cluster_centers_
    labels = kmeans.labels_

    figure(610, figsize=(13,10))
    clf()

    # initialize some arrays:
    etamean_cluster = zeros(n_clusters)
    events_cluster = []
    points_cluster = []
    sample_point_cluster = zeros((n_clusters,2))
    sample_event_cluster = []
    prob_cluster = []
    center_of_mass = zeros((n_clusters,2))
    sample_point_weighted_cluster = zeros((n_clusters,2))
    sample_event_weighted_cluster = []

    for jcluster in range(n_clusters):

        # points (etamean,etamax) in this cluster:
        points_j = cluster_data[labels == jcluster]
        points_cluster.append(points_j)

        # indices of events in this cluster:
        events_j = [j for j in range(400) if labels[j]==jcluster]
        events_cluster.append(events_j)

        # etamean for this cluster, to sort them later by size:
        etamean_cluster[jcluster] = points_j[:,0].mean()

        # sum up probabilities for all events in cluster:
        eprobs = [events_prob[j] for j in events_j]
        prob_cluster.append(sum(eprobs))
                
        # find the event that is closest to the center of the cluster
        # this is the sample event that could be used to represent this cluster
        # in doing fine grid runs:
        sample_distance = inf
        for event,point in zip(events_j,points_j):
            distance_from_center = sqrt((point[0] - centers[jcluster][0])**2 + \
                    (point[1] - centers[jcluster][1])**2)
            if distance_from_center < sample_distance:
                sample_distance = distance_from_center
                sample_point = point
                sample_event = event
        sample_point_cluster[jcluster] = sample_point
        sample_event_cluster.append(sample_event)


        # find the event that is closest to the center of mass of the cluster,
        # weighted by the probability of events:
        sample_distance = inf
        
        for i in range(len(eprobs)):
            center_of_mass[jcluster,:] += points_j[i,:]*eprobs[i]
        center_of_mass[jcluster,:] /= sum(eprobs)

        for event,point in zip(events_j,points_j):
            distance_from_center = sqrt((point[0] - center_of_mass[jcluster,0])**2 + \
                    (point[1] - center_of_mass[jcluster,1])**2)
            if distance_from_center < sample_distance:
                sample_distance = distance_from_center
                sample_point = point
                sample_event = event
        sample_point_weighted_cluster[jcluster] = sample_point
        sample_event_weighted_cluster.append(sample_event)
            

    ####
    # Use the weighted clusters:
    sample_point_cluster = sample_point_weighted_cluster
    sample_event_cluster = sample_event_weighted_cluster
    ####
        
    jsort = etamean_cluster.argsort()   # sort order of index for increasing mean


    cluster_colors = 20*['b','r','g','m','y']

    for jcluster in range(n_clusters-1,-1,-1):
        j = jsort[jcluster]
        plot(points_cluster[j][:,0], points_cluster[j][:,1], \
            'o', color=cluster_colors[jcluster], label = 'cluster %s' % j, \
            markersize=4)
        plot(sample_point_cluster[j,0], sample_point_cluster[j,1], \
            'o',  color=cluster_colors[jcluster], markersize=8)
        plot(centers[j,0], centers[j,1], '+', color=cluster_colors[jcluster], \
            markersize=15, linewidth=4)

        if 0:
            # weighted:
            plot(sample_point_weighted_cluster[j,0], sample_point_weighted_cluster[j,1], \
                'd',  color=cluster_colors[jcluster], markersize=7)
            plot(center_of_mass[j,0], center_of_mass[j,1], 'x', \
                color=cluster_colors[jcluster], markersize=15)

    #legend(loc='lower right', fontsize=10)
    xlabel(v1label,fontsize=20)
    ylabel(v2label,fontsize=20)
    title(plot_title, fontsize=20)

    if save_png:
        fname = 'clusters_%s.png' % scenario_name
        savefig(fname, bbox_inches='tight')
        print 'Created ', fname

    # Plot probabilities of each event and of clustered events:

    figure(200, figsize=(9,6))
    clf()
    for j,Mw in enumerate(Mwlist):
        semilogy(events_etamean[j*100:(j+1)*100], events_prob[j*100:(j+1)*100], \
            'o', color=Mw_colors[Mw], label="Mw=%g" % Mw, markersize=4)
        title('Probabilities of 400 events')
        xlabel('etamean on coarse grid')
        ylabel('probability of event')
        ylim(1e-8,1)
    legend(loc='lower left')

    if save_png:
        fname = 'probabilities_all_runs.png'
        savefig(fname, bbox_inches='tight')
        print 'Created ', fname

    figure(201, figsize=(9,6))
    clf()
    for jcluster in range(n_clusters-1,-1,-1):
        j = jsort[jcluster]
        indices = events_cluster[j]
        etamean = [events_etamean[i] for i in indices]
        eprobs = [events_prob[i] for i in indices]
        semilogy(etamean, eprobs, 'o', color=cluster_colors[jcluster], \
                label = 'cluster %s' % j, markersize=4)
        jsample = sample_event_cluster[j]
        semilogy(events_etamean[jsample], events_prob[jsample], 'o', \
                color=cluster_colors[jcluster], markersize=7)
        semilogy(events_etamean[jsample], prob_cluster[j], 's', \
                color=cluster_colors[jcluster], markersize=6)
        
    #legend(loc='lower left')
    #title('Probabilities of 400 events and clusters')
    title('k-mean clustering of events based on %s and %s' \
            % (v1label,v2label))
    xlabel('etamean on coarse grid')
    ylabel('probability of event')
    ylim(1e-8,1)

    if save_png:
        fname = 'probabilities_%s.png' % scenario_name
        savefig(fname, bbox_inches='tight')
        print 'Created ', fname

    print "\nEvents to use and clustered probability:"
    for jcluster in range(n_clusters):
        j = jsort[jcluster]
        indices = events_cluster[j]
        jsample = sample_event_cluster[j]
        print "Event %s: %s, prob = %7.5f, (%s events in cluster)" \
            % (str(jsample).rjust(3), events_description[jsample], \
            prob_cluster[j], str(len(indices)).rjust(3))

    if make_files:

        # Create output files with clustered events:

        fname = 'cluster_%s.txt' % scenario_name
        cluster_file = open(fname, 'w')
        cluster_file.write( "Events to use and clustered probability:\n")
        for jcluster in range(n_clusters):
            j = jsort[jcluster]
            indices = events_cluster[j]
            jsample = sample_event_cluster[j]
            cluster_file.write( "Event %s: %s, prob = %7.5f (%s events in cluster)\n" \
                % (str(jsample).rjust(3),  events_description[jsample], \
                  prob_cluster[j], str(len(indices)).rjust(3)))

        cluster_file.write( "\nEvents in each cluster:\n")
        for jcluster in range(n_clusters):
            j = jsort[jcluster]
            indices = events_cluster[j]
            jsample = sample_event_cluster[j]
            cluster_file.write( "\nCluster %s is represented by:\n" % j)
            cluster_file.write( "   Event %s: %s, cluster prob = %7.5f\n" \
                % (str(jsample).rjust(3),  events_description[jsample], \
                prob_cluster[j]))
            cluster_file.write( "Clustered events:\n")
            for i in indices:
                cluster_file.write( "   Event %s: %s, prob = %7.5f\n" \
                    % (str(i).rjust(3), events_description[i], events_prob[i]))
        cluster_file.close()
        print "\nCreated ", fname

        # Create scenario file with columns   Mw, runno, probability
        # Order by Mw and then runno

        fname = 'scenario_%s.txt' % scenario_name
        scenario_file = open(fname, 'w')
        jsort2 = argsort(sample_event_cluster)
        for jcluster in range(n_clusters):
            j = jsort2[jcluster]
            jsample = sample_event_cluster[j]
            scenario_file.write( "%3.1f  %3i  %10.8f\n" \
                % (events_Mw[jsample], events_runno[jsample], prob_cluster[j]))
        scenario_file.close()
        print "\nCreated ", fname

        # Create files for Loyce's code:

        fname = 'clusters_for_ptha_%s.txt' % scenario_name
        cluster_file = open(fname, 'w')

        for jcluster in range(n_clusters):
            j = jsort[jcluster]
            indices = events_cluster[j]
            jsample = sample_event_cluster[j]

            # cluster_no      number_in_cluster:
            cluster_file.write( "%3i   %3i\n" % (j,len(indices)))

            # centroid_no     centroid_mag*10
            cluster_file.write( "%s\n" % events_description2[jsample])

            for i in indices:
                cluster_file.write( "%s\n" % events_description2[i])
        cluster_file.close()
        print "\nCreated ", fname


#show()
