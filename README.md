# fema_source_filter

Code for project on source filtering funded by FEMA Region IX.

**Final report:** To appear.

## To reproduce results in report:

Skip to Step 3 if you want to used archived GeoClaw results rather than running all the tsunami simulations.

### Step 1:

Use the code in directory `geoclaw_driver` to set up and run all realizations.  

**Note:** 400 tsunami simulations are performed with a fine grid near Crescent City, and then done again with a coarser grid, so running this takes several days of CPU time.

See `geoclaw_driver/README.md` for instructions.

### Step 2:

Use the code in directory `geoclaw_postprocess` to post process the GeoClaw results.  This produces some relatively small `.npy` files with the fgmax results and some other small files that are needed for the steps below.

See `geoclaw_postprocess/README.md` for instructions.

### Step 3:

If you want to download the results for Steps 1 and 2 as computed for this report, you can do:

    source rsync_all_runs.sh

Or, if you ran the GeoClaw simulations on one computer (e.g. in the cloud) and now want to do the analysis and plotting steps below on a different computer (e.g. your laptop), then you can modify the file `rsync_all_runs.sh` to point to the computer where the runs were performed and then do:

    source rsync_all_runs.sh
 
If you rsync results via this step, then they will go into a new directory `geoclaw_output_via_rsync/all_runs_npy_files`

To continue with the steps below, you should move this directory to `geoclaw_output/all_runs_npy_files`
making sure that you don't overwrite your own computed results in the process!

### Step 4:

Use the code in directory `geoclaw_analysis` to plot the results from the GeoClaw runs, to produce scatter plots of various quantities, and to compute clusters of the realizations that are used in the next step.

See `geoclaw_analysis/README.md` for instructions.
 
