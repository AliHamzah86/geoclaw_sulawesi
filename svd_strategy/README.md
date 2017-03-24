# SVD strategy #

The script `svd_coarse_fine.py` uses all coarse grid runs and a few fine grid
runs to construct approximate p-contours and zeta-contours using the
SVD strategy, detailed in the report.

## Creating approximate hazard maps
```
python svd_coarse_fine
```
reads in fixed-grid inundation values stored in the Geoclaw runs lying in
`../geoclaw_output` and performs the approximation.

## Resulting hazard maps and error analysis

The output will be stored in a subdirectory `_output`. This includes
- hazard maps. For example, `zetacontour_approx_0.8.png`
  is the approximate hazard map the zeta contour with
  p=0.8.
- error of the approximation for each hazard map,
  `zetacontour_error_0.8.png` is the contour plot of the
  difference between the true and the approximate map
- coefficient plots `cvf_Usim.png`, `cvf_Ysim.png` illustrating
  the similarity of the coarse and the fine grid SVDs.
- coefficient plots of the transform matrices `Y_transform_8.8.png`
  for each magnitude.
- saved numpy array of all generated maps `hazardcurves_probs.npy`
- summary of errors `output_error.txt`
- pickled fixed-grid data `data.pkl`
