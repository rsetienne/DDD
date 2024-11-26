test_that("optimizer works", {
  brts <- 1:10
  initparsopt <- c(0.3,0.05,12)
  idparsopt <- 1:3
  res_simplex <- dd_ML(brts = brts,
        initparsopt = initparsopt,
        idparsopt = idparsopt,
        optimmethod = 'simplex')
  res_subplex <- dd_ML(brts = brts,
        initparsopt = initparsopt,
        idparsopt = idparsopt,
        optimmethod = 'subplex')
  #res_nloptr <- dd_ML(brts = brts,
  #      initparsopt = initparsopt,
  #      idparsopt = idparsopt,
  #      optimmethod = 'NLOPT_LN_SBPLX')
  testthat::expect_equal(res_simplex$fvalues,res_subplex$fvalues,tolerance = 0.0001)
  #testthat::expect_equal(res_nloptr$fvalues,res_subplex$fvalues,tolerance = 0.0001)
})

  