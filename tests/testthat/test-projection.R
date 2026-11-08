rel_diff = function(a, b){ max(abs(unlist(a) - unlist(b))) / max(abs(unlist(b))) }

test_that("the path transition matrix gives continuous splines", {
  set.seed(7)
  for(V_Smooth in c(FALSE, TRUE)){
    ts = sort(c(0, stats::runif(6, 0, 5), 5))
    Tm = constructUnboundedCubicSplineTransitionMatrix_Path(ts, V_Smooth)
    expect_equal(ncol(Tm), if(V_Smooth) 4 + 6 else 4 + 2*6)
    expect_equal(Tm[1:4, 1:4], diag(4))

    coef = c(Tm %*% stats::rnorm(ncol(Tm)))
    for(j in 2:7){
      prev = coef[4*(j-2) + 1:4]
      cur = coef[4*(j-1) + 1:4]
      h = ts[j] - ts[j-1]
      expect_equal(cur[4], sum(prev * c(h^3, h^2, h, 1)))          # position
      expect_equal(cur[3], sum(prev * c(3*h^2, 2*h, 1, 0)))        # velocity
      if(V_Smooth) expect_equal(cur[2], sum(prev * c(3*h, 1, 0, 0)))  # acceleration
    }
  }
})

test_that("the path transition matrix handles a single region", {
  expect_equal(constructUnboundedCubicSplineTransitionMatrix_Path(c(0, 1), FALSE), diag(4))
  expect_equal(constructUnboundedCubicSplineTransitionMatrix_Path(c(0, 1), TRUE), diag(4))
})

test_that("block-diagonal helpers match the full products", {
  set.seed(8)
  blocks = list(matrix(stats::rnorm(16), 4), matrix(stats::rnorm(16), 4), matrix(stats::rnorm(16), 4))
  M = matrix(stats::rnorm(12*5), 12, 5)
  N = matrix(stats::rnorm(5*12), 5, 12)
  full = VFmixing:::blockDiagDense(blocks)
  expect_equal(VFmixing:::blockDiagTimes(blocks, M), full %*% M)
  expect_equal(VFmixing:::timesBlockDiag(N, blocks), N %*% full)
})

test_that("Projection reproduces the saved example path", {
  # ExEstPath_Projection was produced by the original implementation.
  proj = quietly(getEstimatedPaths_Projection(ExData, 100, V_Smooth = FALSE, Model_Selection_Type = "BIC", undersmooth = FALSE))
  expect_equal(proj$TimeSplits, ExEstPath_Projection$TimeSplits)
  expect_lt(rel_diff(proj$Path, ExEstPath_Projection$Path), 1e-6)
  expect_lt(rel_diff(proj$Cov, ExEstPath_Projection$Cov), 1e-6)
  expect_lt(rel_diff(proj$Free_Parameter_Cov, ExEstPath_Projection$Free_Parameter_Cov), 1e-6)
})
