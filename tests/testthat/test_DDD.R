context("test_DDD")

test_that("DDD works", {
  expect_equal_x64 <- function(object, expected, ...)
  {
    if(Sys.getenv("R_ARCH") == "/x64")
    {
       result <- testthat::expect_equal(object = object, expected = expected,...)
    } else
    {
      result <- NULL
    }
    return(invisible(result))
  }
  
  pars1 = c(0.8,0.1,40)
  pars2 = c(100,1,1,0,0,2)
  brts = 1:30
  missnumspec = 0
  methode = 'lsoda'
  
  r0 <- DDD:::dd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,methode = methode)
  r1 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_rhs',methode = methode)
  r2 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_rhs_FORTRAN',methode = methode)
  r3 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,methode = 'analytical')

  testthat::expect_equal(r0,r2,tolerance = .00001)
  testthat::expect_equal(r1,r2,tolerance = .00001)
  testthat::expect_equal(r1,r3,tolerance = .00001)
  testthat::expect_equal(-39.36047,r1,tolerance = .000001)
  
  r4 <- dd_SR_loglik(pars1 = c(0.2,0.1,50,0.2,0.1,70,5), pars2 = c(100,1,1,1,0,2), brts = 1:10, missnumspec = 0)
  testthat::expect_equal(-27.37304,r4,tolerance = .000001)

  brts = 1:5
  pars2 = c(100,1,3,0,0,2)
  r5 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_bw_rhs',methode = methode)
  r6 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_bw_rhs_FORTRAN',methode = methode)
  r7 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,methode = 'analytical')
  
  testthat::expect_equal(r5,r6,tolerance = .00001)
  testthat::expect_equal(r5,r7,tolerance = .01)
  expect_equal_x64(-8.582413,r7,tolerance = .00001)

  pars1 = c(0.2,0.05,1000000)
  pars2 = c(1000,1,1,0,0,2)
  brts = 1:10
  r8 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_rhs_FORTRAN',methode = methode)
  r9 <- DDD:::dd_loglik_test(pars1 = c(pars1[1:2],Inf),pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_rhs_FORTRAN',methode = methode)
  expect_equal_x64(r8,r9,tolerance = .00001)
  
  pars1 <- c(0.2,0.05,15)
  pars2 <- c(4,1,0,0,2)
  brts <- 1:10
  r10 <- DDD::bd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = 0,methode = 'lsoda')
  expect_equal_x64(r10,-11.6076802136507471,tolerance = 1E-16)

  r11 <- DDD::bd_loglik(pars1 = c(0.4,0.1,20),pars2 = c(4,0,1,0,2), brts = 1:10, missnumspec = 0)
  expect_equal_x64(r11,-27.7337684064852610,tolerance = 1E-16)
  
  pars <- c(0.4, 0.1, 0.6, 0.2)
  brts <- list(10:6, 3:2)
  #sls_test <- sls::loglik_sls_p(
  #  pars = pars,
  #  brts = brts,
  #  cond = 0,
  #  n_max = 1e3
  #);
  high_k <- 1e7
  pars1 <- c(pars[1], pars[2], high_k, pars[3], pars[4], high_k, brts[[2]][1])
  pars2 <- c(200,1,0,brts[[2]][1],0,2,1.5)
  brtsM <- brts[[1]]
  brtsS <- brts[[2]][-1]
  ddd_test <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test,-24.4172040142037936,tolerance = .000001)
  
  pars <- c(0.35, 0.15, 0.7, 0.25)
  brts <- list(c(10:5,2), c(4,3,1))
  # sls_test <- sls::loglik_sls_p(
  #   pars = pars,
  #   brts = brts,
  #   cond = 0,
  #   n_max = 1e3
  # );
  t_d <- brts[[2]][1]
  tsplit <- min(abs(brts[[1]][abs(brts[[1]]) > t_d]))
  high_k <- 1e7
  pars1 <- c(pars[1], pars[2], high_k, pars[3], pars[4], high_k, t_d)
  pars2 <- c(200,1,0,brts[[2]][1],0,2,1.5)
  brtsM <- brts[[1]]
  brtsS <- brts[[2]][-1]
  ddd_test <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'analytical'
  );
  testthat::expect_equal(ddd_test,-29.5449889636434051,tolerance = .000001)
  
  brts <- 1:10
  pars1 <- c(0.8,0,40)
  pars2 <- c(2,1,1,1,1)
  missnumspec <- 0
  result1 <- DDD::bd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec)
  result2 <- DDD::dd_loglik(pars1 = pars1,pars2 = c(min(1000,10 * (length(brts) + missnumspec)),1,pars2[2:5]),brts = brts,missnumspec = 0)
  testthat::expect_equal(result1,result2)
})


