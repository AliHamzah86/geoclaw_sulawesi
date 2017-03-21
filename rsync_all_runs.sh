
# Retrieve all_runs_npy_files if you performed the GeoClaw runs on a
# different computer than where you plan to do the analysis.  

# Change the username, computer address, and directory to point to
# geoclaw_output/all_runs_npy_files.

# If you don't change them, this will pull down results from the 
# runs done at UW and used to produce the final report:

# After running this, move geoclaw_output_via_rsync/all_runs_npy_files
# to geoclaw_output/all_runs_npy_files  for further analysis scripts to work

mkdir -p geoclaw_output_via_rsync

rsync -avz \
    ptha@homer.u.washington.edu:public_html/fema2016/all_runs_npy_files \
    geoclaw_output_via_rsync/


