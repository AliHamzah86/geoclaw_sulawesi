
# GeoClaw post-processing scripts

These scripts should be run on the same computer where the GeoClaw runs were performed.

make_npy_files_Mw.py 
    
    Creates geoclaw_output/all_runs_npy_files 
    
    This produces fgmax output, rewritten as compact npy files.
    Could also make gauge npy files, currently not used.
    

compute_dtopo_proxies.py    

    Reads each dtopo file e.g. from geoclaw_output/coarse_8.6/run_*
    and produces table for all runs at this magnitude in file:
    geoclaw_output/all_runs_npy_files/dtopo_proxies_8.6.txt

make_B0.py
    Creates files 
        geoclaw_output/fine_xyB0_fg1.npy  etc.
    Uses fgmax output from 
        geoclaw_output/fine_B0
        geoclaw_output/coarse_B0
        
    Note: Makefile should be set to use 
        $(GEOLIB)/fgmax_interpolate0.f90
    from Version 5.3.1 to get the essentially the same B0 as used in the report.
    For the work done in the report, B0 was calculated from a run that had 
    essentially zero subsidence at Crescent City, but it varied slightly over 
    the fgmax domain.  
    
