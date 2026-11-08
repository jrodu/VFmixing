rel_diff = function(a, b){ max(abs(a - b)) / max(abs(b)) }

test_that("quadrature integrals match the adaptive integrals", {
  P = ExEstPath_FullFours
  ts = c(0, 1.7, 5)   # trajectory nodes not aligned with the path nodes
  gl = VFmixing:::gaussLegendre(20)

  adaptive_traj = constructIntegralList_Traj(P, baseVectorFields, list(), ts, 1, 2)
  adaptive_pathtraj = constructIntegralList_PathTraj(P, baseVectorFields, list(), 1, P$TimeSplits, ts, 2, 2)

  for(r in 1:2){
    cache = VFmixing:::oneShotRegionCache(ts[r], ts[r + 1], P, baseVectorFields, 2, 2, gl)
    expect_lt(rel_diff(cache$G, adaptive_traj[[r]]), 1e-8)
    expect_lt(rel_diff(t(cache$S), adaptive_pathtraj[[r]]), 1e-8)
  }
})

test_that("Gauss-Legendre rule integrates polynomials exactly", {
  gl = VFmixing:::gaussLegendre(5)
  for(p in 0:9){
    expect_equal(sum(gl$w * gl$x^p), if(p %% 2 == 0) 2/(p + 1) else 0, tolerance = 1e-12)
  }
})

test_that("OneShot reproduces the saved example trajectory", {
  # A regression check against ExEstTraj_OneShot. The comparison with the original
  # adaptive integrals is the "quadrature integrals match" test above.
  fast = quietly(getEstimatedTrajectories_OneShot(ExEstPath_FullFours, nsplits = 100, baseVectorFields = baseVectorFields,
                                                  undersmooth = FALSE))

  expect_equal(fast$TimeSplits, ExEstTraj_OneShot$TimeSplits)
  expect_lt(rel_diff(unlist(fast$Traj), unlist(ExEstTraj_OneShot$Traj)), 1e-6)
  expect_lt(rel_diff(unlist(fast$Cov), unlist(ExEstTraj_OneShot$Cov)), 1e-6)
  expect_lt(rel_diff(fast$Error, ExEstTraj_OneShot$Error), 1e-6)
  expect_setequal(names(fast), names(ExEstTraj_OneShot))
})

test_that("patience stops early but keeps the best model when it comes first", {
  full = quietly(getEstimatedTrajectories_OneShot(ExEstPath_FullFours, nsplits = 30, baseVectorFields = baseVectorFields))
  early = quietly(getEstimatedTrajectories_OneShot(ExEstPath_FullFours, nsplits = 30, baseVectorFields = baseVectorFields, patience = 5))
  expect_equal(early$TimeSplits, full$TimeSplits)
  expect_equal(early$Traj, full$Traj)
})

test_that("OneShot validates quad_nodes and patience", {
  expect_error(getEstimatedTrajectories_OneShot(ExEstPath_FullFours, 5, baseVectorFields, quad_nodes = 0), "quad_nodes")
  expect_error(getEstimatedTrajectories_OneShot(ExEstPath_FullFours, 5, baseVectorFields, quad_nodes = 2.5), "quad_nodes")
  expect_error(getEstimatedTrajectories_OneShot(ExEstPath_FullFours, 5, baseVectorFields, patience = 0), "patience")
  expect_error(getEstimatedTrajectories_OneShot(ExEstPath_FullFours, 5, baseVectorFields, Model_Selection_Type = "nope"), "must be one of")
})

test_that("linear-algebra helpers match the originals", {
  ts = c(0, 1, 2.5, 3, 5)
  expect_equal(VFmixing:::trajTransitionMatrix(ts, 2, FALSE), constructUnboundedCubicSplineTransitionMatrix_Traj(ts, 2, FALSE))
  expect_equal(VFmixing:::trajTransitionMatrix(ts, 3, TRUE), constructUnboundedCubicSplineTransitionMatrix_Traj(ts, 3, TRUE))

  set.seed(1)
  X = matrix(stats::rnorm(40), 10, 4)
  B = matrix(stats::rnorm(8), 4, 2)
  pd = crossprod(X)
  expect_equal(VFmixing:::psdSolve(pd, B), pracma::pinv(pd) %*% B)

  singular = crossprod(cbind(X[, 1:3], X[, 1] + X[, 2]))   # rank 3
  expect_equal(VFmixing:::psdSolve(singular, B), pracma::pinv(singular) %*% B)
})
