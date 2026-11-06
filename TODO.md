# VFmixing TODO

Notes on things to work on later. Check items off (`- [x]`) as they're done.

## Open

_Nothing open right now. Add new items here._

## Done

### Bugs

- [x] **`getBestSplitsTraj` `"BICPos"` can't run.** Removed: it needed a velocity
  covariance matrix that no trajectory method produces. Asking for it now stops
  with a message pointing to `"BIC"` or `"VarPos"`.
- [x] **`getBestSplitsTraj` `"BIC"` and `"AIC"` assume `Error` is a matrix.**
  Both now use one error summary that handles ByFours (matrix) and OneShot
  (per-region vector) trajectories. `"BICOS"` is kept as an alias of `"BIC"`.
  `"VarPos"` stops clearly if the trajectory has no covariance.
- [x] **`getBestSplitsPath` `"BICVel"` is wrong in two places.** Fixed the
  parameter count and the indexing, plus a third bug: the variance-jump check
  read a column (`BIC_Vel[,2]`) instead of the variance row (`BIC_Vel[2,]`).
  If the variance never jumps, every model is considered.
- [x] **`coverageRateTraj` returns one pooled rate.** Now returns one rate per model.

### Consistency

- [x] **`n_splits` means different things.** Every method now stores the number
  of regions. The "Best Number of Splits" message is now "Best number of regions".
- [x] **OneShot uses one `V_Smooth` for both path and trajectory.** The path's
  smoothness is now detected from its free parameters. ByFours trajectories do
  the same when `Path_V_Smooth` is not given.
- [x] **Default arguments differ between path methods.** Projection now defaults
  to `V_Smooth = FALSE`, `"BIC"`, like FullRegression and FullFours.
- [x] **`drifterVisualizationTool`: `requireNamespace()` only checks its first
  argument.** Now checks every package and lists the missing ones.
- [x] **2-D assumptions.** `visualizeData` takes `x_var`/`y_var`;
  `animatePathTraj` stops clearly for non-2-D paths; `getPropagationList` works
  in any dimension (2-D keeps the Lon/Lat names).
- [x] **`animatePathTraj` uses the deprecated ggplot2 `size` argument.** Now
  `linewidth`; DESCRIPTION requires ggplot2 >= 3.4.0.
- [x] **`RColorBrewer` is only in Suggests.** Moved to Imports.

### Features

- [x] **`Model_Selection_Type` for FullFours.** `"BIC"`/`"AIC"` select among the
  ByFours candidates directly (fast). `"BICVel"`, `"VarPos"`, `"VarVel"` and
  `"VarPosVel"` refit the full regression on every candidate and select among
  those. ByFours itself now also accepts `"AIC"`.

### Data and docs

- [x] **Regenerate `ExEstTraj_ByFours` and `ExPropagationList`.** Also regenerated
  `ExEstPath_ByFours` for the new `n_splits`. `ExEstTraj_ByFours` went from 42 to
  3 regions under the current BIC.
- [x] **Add a `tests/testthat` suite.**
- [x] **Slow examples.** The OneShot and Projection examples are in `\donttest{}`.
