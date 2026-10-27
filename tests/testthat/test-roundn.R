# roundn rounds half up. On develop (fae24c7) it uses trunc(fac * x + 0.5) / fac, which is
# wrong for negative numbers (roundn(-2) = -1, roundn(-0.5, 2) = -0.49, roundn(-10.2) = -9).

test_that("roundn works for positive numbers", {
  expect_equal(DDD::roundn(2), 2)
  expect_equal(DDD::roundn(2.4), 2)
  expect_equal(DDD::roundn(2.5), 3)
  expect_equal(DDD::roundn(0.5), 1)
  expect_equal(DDD::roundn(1.2345, 3), 1.235)
})

test_that("roundn works for negative numbers", {
  expect_equal(DDD::roundn(-2), -2)
  expect_equal(DDD::roundn(-0.5, 2), -0.5)
  expect_equal(DDD::roundn(-10.2), -10)
  expect_equal(DDD::roundn(-2.4), -2)
  expect_equal(DDD::roundn(-1.2345, 3), -1.234)
})
