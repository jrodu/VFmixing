# VFmixing TODO

Notes on things to work on later. Check items off (`- [x]`) as they're done.

## Bugs

- [ ] **`getBestSplitsTraj` `"BICPos"` can't run.** It reads `cur_traj$EstVelCovMat`,
  which no trajectory estimator produces, and indexes `Error` as a matrix.
  Decide what BICPos should compute, or remove it.
- [ ] **`getBestSplitsTraj` `"BIC"` and `"AIC"` assume `Error` is a matrix.**
  OneShot's `Error` is a vector of per-region integrals, so `Error[,-1]` fails.
  `"AIC"` also takes the median of `Error[,1]`, which is the time column.
- [ ] **`getBestSplitsPath` `"BICVel"` is wrong in two places.**
  It counts parameters with `length(cur_path)` (the number of list fields)
  instead of `length(cur_path$Path)`, and indexes the 2 x K result matrix
  linearly: `BIC_Vel[1:vel_cutoff]` should be `BIC_Vel[1, 1:vel_cutoff]`.
  `vel_cutoff` is also `NA` if the variance never jumps.
- [ ] **`coverageRateTraj` returns one pooled rate.** The docs say one rate per
  model; the `matrix(...)` wrapper flattens the result before `colMeans`.

## Consistency

- [ ] **`n_splits` means different things.** ByFours stores the number of interior
  nodes; FullRegression, Projection and FullFours store the number of regions.
- [ ] **OneShot uses one `V_Smooth` for both path and trajectory.** If the path
  was fit with different smoothness, the free-parameter dimensions won't match.
- [ ] **Default arguments differ between path methods.** Projection defaults to
  `V_Smooth = TRUE`, `"VarPosVel"`; FullRegression to `FALSE`, `"BIC"`.
- [ ] **`drifterVisualizationTool`: `requireNamespace()` only checks its first
  argument,** so only `leaflet` is actually checked.
- [ ] **2-D assumptions.** `visualizeData` and `animatePathTraj` assume two
  dimensions; `getPropagationList` hard-codes Lon/Lat column names.
- [ ] **`animatePathTraj` uses the deprecated ggplot2 `size` argument** for lines
  (should be `linewidth`).
- [ ] **`RColorBrewer` is only in Suggests** but `visualizeEstTraj` uses it
  without checking it's installed.

## Features

- [ ] **`Model_Selection_Type` for FullFours.** Node selection is currently BIC
  only (hard-coded in ByFours). Options: allow AIC, or refit the full regression
  on several ByFours candidates and compare with the variance-based criteria.

## Data and docs

- [ ] **Regenerate `ExEstTraj_ByFours` and `ExPropagationList`** with the current
  code. `ExPropagationList` likely still comes from the pre-rename trajectory.
- [ ] **Add a `tests/testthat` suite.** There are no tests; `R CMD check` only
  runs the examples.
- [ ] **Slow examples:** `getEstimatedTrajectories_OneShot` (~30 s) and
  `getEstimatedPaths_Projection` (~10 s). Shorten or wrap in `\donttest{}`
  if the package goes to CRAN.
