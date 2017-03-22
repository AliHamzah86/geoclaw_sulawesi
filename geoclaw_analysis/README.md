
# GeoClaw analysis scripts

These scripts assume `geoclaw_output/all_runs_npy_files` exists.

## Compute etamean, etamax, etc. from each run

### `compute_etamax.py`

Produces files such as:
`all_runs_npy_files/coarse_9.0_etamax_fg1.txt`

Containing:

    run     min eta    mean eta    max eta (on shore)  hmax sum
     0       6.885       7.087       7.314              186.911
     1       5.623       5.893       6.137              105.705
     etc.
     
## Plot fgmax from each run and scatter plots comparing things

### `plot_fgmax_transects_all_runs_Mw.py`

Make plots of fgmax grids 1 and 2 (1d transects), putting 50 runs on the plot for each Mw / resolution.

Creates files e.g. `png_files/Mw86_fg1_50runs_coarse.png`
     
### `plot_fgmax_grid_all_runs_Mw.py`

Makes plots of fgmax grid 3 (2d grid) for each run individually.

Creates files e.g. `png_files/coarse_9.0/run_1/run1_Mw9.0_fg3.png`

### `plot_coarse_fine_scatter.py`

 Plots fine vs. coarse scatter plots of etamean or other qoi
 Plots etamean vs. etamax for each resolution.  
 Puts these figures in `png_files`

     Created  scatter_etamean_fg1.png
     Created  scatter_etamean_fg2.png
     Created  scatter_etamean_fg3.png

### `plot_scatter_for_report.py`

Creates scatter plot used in report.  Puts these figures in `png_files`
      
      Created  scatter_etamean_fg3.png
      Created  scatter_etamean_fg3_zoom.png
      Created  scatter_etamax_fg3.png
      Created  scatter_etamax_fg3_zoom.png
      Created  scatter_volume_fg3.png
      Created  scatter_volume_fg3_zoom.png
      Created  scatter_etamean_etamax_coarse_fg3.png
      Created  scatter_etamean_etamax_fine_fg3.png

## Create clusters

### `kmeans_cluster.py`

k-means clustering based on cluster_vars, e.g.
dtopo_proxies or coarse_mod_eta

Requires `../scenario_prb_wgts.txt`

## Scripts to make html and latex files

### `collect_all_realizations.py`

Make html file of all realizations slip, dtopo, hmax together

