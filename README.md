
<!-- README.md is generated from README.Rmd. Please edit that file -->

# VFmixing

<!-- badges: start -->
<!-- badges: end -->

The goal of VFmixing is to provide methodology specifically designed to
mix vector fields given noisy positional data.

## Installation

You can install the development version of VFmixing from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jrodu/VFmixing")
```

## Example

This example shows how to use the main set of functions. Let’s look at
the noisy positional data from a particle moving through some unknown
vector field.

``` r
library(VFmixing)

Data <- VFmixing::ExData
visualizeData(Data)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

The key assumption for this method is that the true unknown vector field
is a linear combination of a set of simpler known vector fields at every
point in time. For this example, the set of mixing vector fields is a
rotational field and an expanding field. To use this set of vector
fields in the mixing process we construct the following function:

``` r
baseVectorFields
#> function (t, curPos) 
#> {
#>     f1 = c(curPos[2], -1 * curPos[1])/sqrt(sum(c(curPos[1], curPos[2])^2))
#>     f2 = c(curPos[1], curPos[2])/sqrt(sum(c(curPos[1], curPos[2])^2))
#>     matrix(c(f1, f2), nrow = 2, byrow = F)
#> }
#> <bytecode: 0x000001b589766508>
#> <environment: namespace:VFmixing>
```

Now that the `baseVectorFields` function has been defined, we can begin
the mixing process. This is broken into two steps: Path and Trajectory.

### Step 1: Estimate the path

The first step estimates the path the particle took through the vector
field. We do this to get better estimates of both position and velocity,
which would otherwise be compromised extremely quickly with only a small
degree of error. The path is estimated with a cubic spline using
`getEstimatedPaths_FullFours()`.

The “FullFours” method works in two stages:

1.  The fast “ByFours” algorithm searches for the spline nodes. Each
    time a region is split, only the four neighboring segments are
    refit, so the search is extremely fast. By default, the number of
    nodes is chosen by BIC.
2.  A single full regression is then fit on the selected nodes. This
    gives stable estimates and an estimable covariance structure for the
    path.

The result is the speed of “ByFours” together with the covariance
structure needed for confidence intervals, random path draws, and the
“One-Shot” trajectory method below.

The criterion used to choose the nodes can be changed with
`Model_Selection_Type`. `"BIC"` (the default) and `"AIC"` are applied
directly to the “ByFours” candidates, so they stay fast. The
covariance-based criteria (`"BICVel"`, `"VarPos"`, `"VarVel"` and
`"VarPosVel"`) need a full regression on every candidate, so they are
slower. Every path method also takes a `patience` argument, which stops
the search once that many candidates in a row fail to improve the
criterion (see Step 2 for an example). It makes the biggest difference
with the slower criteria and with `getEstimatedPaths_FullRegression()`.

``` r
PathTime <- system.time(
  EstPath <- getEstimatedPaths_FullFours(Data = Data, nsplits = 100)
)
```

The run time, in seconds:

``` r
PathTime
#>    user  system elapsed 
#>    0.18    0.08    0.32
```

Because the “FullFours” path has a covariance structure, the plot can
show its uncertainty. By default these are 95% prediction intervals for
a new observation: the uncertainty in the fitted path plus the
observation noise (`P_Int = FALSE` shows a confidence interval for the
path itself instead). When both axes are positions, the shaded tube is
the joint 95% region for both coordinates, swept along the path:

``` r
visualizeEstPath(Data = Data, EstPath = EstPath, t_grid_size = 0.01,
                 x_var = 'X1', y_var = 'X2')
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

With time on one axis, the interval is drawn as a band:

``` r
visualizeEstPath(Data = Data, EstPath = EstPath, t_grid_size = 0.01,
                 x_var = 't', y_var = 'X1')
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

We can also check how well this uncertainty is calibrated.
`coverageRatePath()` uses the same 95% prediction intervals at each data
time and reports, for each dimension, the share of data points that fall
inside them.

``` r
PathCovRate = coverageRatePath(EstPath = EstPath, Data = Data, CL = 0.95)
PathCovRate
#> [1] 0.9441118 0.9620758
```

The coverage rates for X1 and X2 are 94.4% and 96.2%, close to the
nominal 95%.

This estimated path is also included in the package as
`ExEstPath_FullFours`.

<details>
<summary>
Other path methods
</summary>

The package also provides the methods that “FullFours” is built from.
They are mainly useful for comparison:

- `getEstimatedPaths_ByFours()`: the node search on its own. Very fast,
  but without a covariance structure.
- `getEstimatedPaths_FullRegression()`: refits the full spline at every
  split. Has a covariance structure, but is much slower.
- `getEstimatedPaths_Projection()`: fits each region separately and
  projects onto the smooth spline space. Has a covariance structure, but
  estimates can be unstable for very complex data.

</details>

### Step 2: Estimate the trajectories

With the path of the particle estimated, we can move to the main event:
the trajectory step. We want to estimate the weight functions for each
mixing vector field. We call the weight functions the Trajectory of each
field. These are also estimated with cubic splines.

Because the “FullFours” path has a covariance structure, we can use the
“One-Shot” method. It uses the path covariance structure to transform
the path coefficients into the trajectory coefficients, which yields a
theoretical covariance structure for the trajectories under the
assumption that the mixing vector fields are relatively constant in a
small window defined by the path error. A faster alternative without a
covariance structure is `getEstimatedTrajectories_ByFours()`.

``` r
TrajTime <- system.time(
  EstTraj <- getEstimatedTrajectories_OneShot(EstPath = EstPath, nsplits = 100,
                                              baseVectorFields = baseVectorFields)
)
```

The run time, in seconds:

``` r
TrajTime
#>    user  system elapsed 
#>   11.74    0.11   12.17
```

By default, every candidate up to `nsplits` is fit before the best one
is chosen. The `patience` argument stops the search once that many
candidates in a row have failed to improve the selection criterion. This
is much faster when the best model has few regions, but it can miss a
better model that would only appear after a long run without
improvement, so it is off by default.

``` r
PatienceTime <- system.time(
  EstTrajPatience <- getEstimatedTrajectories_OneShot(EstPath = EstPath, nsplits = 100,
                                                      baseVectorFields = baseVectorFields,
                                                      patience = 10)
)
```

``` r
PatienceTime
#>    user  system elapsed 
#>    0.08    0.00    0.08

# Same trajectory as the full search?
isTRUE(all.equal(EstTrajPatience$Traj, EstTraj$Traj))
#> [1] TRUE
```

``` r
TrueTraj <- VFmixing::ExTraj
visualizeEstTraj(EstTraj = EstTraj, TrueTraj = TrueTraj$Traj,
                 TrueTrajTimeSplits = TrueTraj$TimeSplits, EstPath = EstPath,
                 baseVectorFields = baseVectorFields, t_grid_size = 0.01,
                 CL = 0.95, P_Int = T, Cos_Cutoff = 0.91)
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

In the plot above, we can see our estimate of the trajectory (solid
line) next to the true trajectories used to construct the data (dashed
line). We can evaluate the method’s performance more formally using the
average squared distance between the estimated trajectories and the
truth, as well as the coverage rate of each trajectory’s 95% intervals.

``` r
AvgSqDist = cubicAvgSquareDist(EstCubic = EstTraj$Traj,
                               TrueCubic = TrueTraj$Traj,
                               EstTimeSplits = EstTraj$TimeSplits,
                               TrueTimeSplits = TrueTraj$TimeSplits,
                               lower_bound = min(EstTraj$TimeSplits),
                               upper_bound = max(EstTraj$TimeSplits))
AvgSqDist
#> [1] 0.0002546923

TrajCovRate = coverageRateTraj(EstTraj = EstTraj, TrueTraj = TrueTraj$Traj,
                               TrueTrajTimeSplits = TrueTraj$TimeSplits,
                               CL = 0.95)
TrajCovRate
#> [1] 1 1
```

The average squared distance is extremely small, and the coverage rates
for the two trajectories are 100% and 100%, at or above the nominal 95%.
For this example the intervals are, if anything, conservative.
