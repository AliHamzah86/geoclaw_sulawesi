# Generate GeoClaw output to be used for source filtering #

This directory contains python scripts that will perform 802 geoclaw runs,
each with different earthquake events in the southern portion of the
Cascadia Subduction Zone, ranging from magnitude 8.6 to 9.2.

## Generating random earthquakes

Parameters for the random earthquake events generated for this report are
stored in `scenario_pts.txt` and `scenario_prb_wgts.txt`.
The code that generated these events are demonstrated in the IPython notebook
`Scenario_generation_PTHA.ipynb`.

## Install GeoClaw

See http://www.clawpack.org/installing.html for general instructions to install Clawpack, 
which includes GeoClaw.  Version 5.4.0 was used for the original runs.

## Adjust the Makefile

Set the `FFLAGS` in the `Makefile` to be appropriate for your compiler.
You probably want to use OpenMP parallelism if possible.

## Fetch topography DEMs

In this directory:

```
$ source get_topo.sh
```

## Perform Geoclaw runs

**Note:** This will perform 400 coarse grid and 400 fine grid runs and will take
several days or weeks to run depending on how many cores are being used...

To start the run
```
$ python run_CC_CSZ_South.py
```
The outputs will be organized in `../geoclaw_output` directory.
For example, magnitude 8.8 event number 10 computed with a coarse grid
will be stored in
```
../geoclaw_output/coarse_8.8/run_10
```

## Generate plots for each event
To create plots for each of these runs, do
```
$ python make_plots_all.py
```
