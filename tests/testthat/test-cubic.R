test_that("evaluateCubic and its derivatives match a known polynomial", {
  # f(t) = 2(t-1)^3 - (t-1)^2 + 3(t-1) + 5
  cubic = c(2, -1, 3, 5)
  t = c(0, 1, 2.5)
  s = t - 1

  expect_equal(c(evaluateCubic(t, cubic, 1)), 2*s^3 - s^2 + 3*s + 5)
  expect_equal(c(evaluateCubicVelocity(t, cubic, 1)), 6*s^2 - 2*s + 3)
  expect_equal(c(evaluateCubicAcceleration(t, cubic, 1)), 12*s - 2)
  expect_equal(c(evaluateCubicJerk(t, cubic, 1)), rep(12, 3))
})

test_that("evaluateCubic returns one column per cubic", {
  cubics = matrix(c(1, 1, 1, 1,
                    1, 2, 3, 4), nrow = 2, byrow = TRUE)
  out = evaluateCubic(0:3, cubics, 0)
  expect_equal(dim(out), c(4, 2))
})

test_that("generateC1Smoothness joins cubics with matching value and slope", {
  prev = matrix(c(1, 2, 3, 4), nrow = 1)
  cur = generateC1Smoothness(prev, matrix(c(2, 0.5), nrow = 1), prevStartTime = 0, curStartTime = 1, meetTime = 1)

  expect_equal(c(evaluateCubic(1, cur, 1)), c(evaluateCubic(1, prev, 0)))
  expect_equal(c(evaluateCubicVelocity(1, cur, 1)), c(evaluateCubicVelocity(1, prev, 0)))
})
