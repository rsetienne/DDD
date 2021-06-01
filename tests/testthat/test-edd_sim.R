testthat::test_that("per species rates check works", {
  testthat::expect_error(
    edd_pars_check(
      pars = c(-0.5, 0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "none"
    ),
    "per species rates should be positive"
  )
  
  testthat::expect_error(
    edd_sim(
      pars = c(-0.5, -0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "none"
    ),
    "per species rates should be positive"
  )
  
  testthat::expect_error(
    edd_sim(
      pars = c(0.5, -0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "none"
    ),
    "per species rates should be positive"
  )
})

testthat::test_that("model and parameters match check works", {
  testthat::expect_error(
    edd_sim(
      pars = c(0.5, 0.1, -0.001, -0.001),
      age = 3,
      model = "dsde2",
      metric = "ed",
      offset = "none"
    ),
    "incorrect parameters"
  )
  
  testthat::expect_error(
    edd_sim(
      pars = c(-0.5, -0.1, -0.001, -0.001, 0.001, 0.001),
      age = 3,
      model = "dsce2",
      metric = "ed",
      offset = "none"
    ),
    "this model requires four parameters"
  )
})