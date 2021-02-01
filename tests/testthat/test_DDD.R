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
  
  r4 <- DDD::dd_SR_loglik(pars1 = c(0.2,0.1,50,0.2,0.1,70,5), pars2 = c(100,1,1,1,0,2), brts = 1:10, missnumspec = 0)
  testthat::expect_equal(-27.37304,r4,tolerance = .000001)

  brts = 1:5
  pars2 = c(100,1,3,0,0,2)
  r5 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_bw_rhs',methode = methode)
  r6 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,rhs_func_name = 'dd_loglik_bw_rhs_FORTRAN',methode = methode)
  r7 <- DDD:::dd_loglik_test(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,methode = 'analytical')
  
  testthat::expect_equal(r5,r6,tolerance = .00001)
  testthat::expect_equal(r5,r7,tolerance = .01)
  expect_equal_x64(-8.579058,r7,tolerance = .00001) #was -8.582413 before

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
})

context("test_DDD_KI")

test_that("DDD_KI works",
{
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
  pars2 <- c(500,1,0,brts[[2]][1],0,2,1.5)
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
  testthat::expect_equal(ddd_test,-24.4171970357049624,tolerance = .000001)
  ddd_test2 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'analytical'
  ); 
  testthat::expect_equal(ddd_test,ddd_test2,tolerance = .000001)
  
  low_k <- 20
  pars1[3] <- low_k
  ddd_test <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test,-21.1781625797899231,tolerance = .000001)
  ddd_test2 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'analytical'
  ); 
  testthat::expect_equal(ddd_test,ddd_test2,tolerance = .000001)
  
  ddd_test3 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 3,
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test3,-19.6273910107265408,tolerance = .000001)
  
  ddd_test03 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = c(0,3),
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test03,-21.4981352311200595,tolerance = .000001)
  
  ddd_test12 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = c(1,2),
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test12,-20.7167138427128776,tolerance = .000001)
  
  ddd_test21 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = c(2,1),
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test21,-20.0405933298708874,tolerance = .000001)
  
  ddd_test30 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = c(3,0),
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test30,-19.4834201422017124,tolerance = .000001)
  
  #testthat::expect_equal(ddd_test3,log(exp(ddd_test03) + exp(ddd_test12) + exp(ddd_test21) + exp(ddd_test30)),tolerance = .000001)
  #not identical due to additional combinatorial factors
  
  pars2[3] <- 1
  ddd_test <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test,-20.5299241171281643,tolerance = .000001)
  ddd_test2 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'analytical'
  ); 
  testthat::expect_equal(ddd_test,ddd_test2,tolerance = .000001)
  
  pars2[3] <- 4
  ddd_test <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test,-20.2509115267895616,tolerance = .000001)

  pars2[3] <- 5
  ddd_test <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'ode45'
  ); 
  testthat::expect_equal(ddd_test,-20.1686905596579997,tolerance = .000001)

  cond <- 1
  pars <- c(0.35, 0.15, 0.7, 0.25)
  brts <- list(c(10:5,2), c(4,3,1))
  # sls_test <- sls::loglik_sls_p(
  #   pars = pars,
  #   brts = brts,
  #   cond = cond,
  #   n_max = 1e3
  # );
  t_d <- brts[[2]][1]
  tsplit <- min(abs(brts[[1]][abs(brts[[1]]) > t_d]))
  high_k <- 1e7
  pars1 <- c(pars[1], pars[2], high_k, pars[3], pars[4], high_k, t_d)
  pars2 <- c(200,1,cond,brts[[2]][1],0,2,1.5)
  brtsM <- brts[[1]]
  brtsS <- brts[[2]][-1]
  ddd_test1 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'analytical'
  );
  testthat::expect_equal(ddd_test,-28.5415506633517460,tolerance = .000001)
  cond <- 0
  pars2 <- c(200,1,cond,brts[[2]][1],0,2,1.5)
  ddd_test0 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'analytical'
  )
  testthat::expect_equal(ddd_test0,-29.5449838542774543,tolerance = .000001)
  ddd_test2 <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = 'ode45'
  )
  testthat::expect_equal(ddd_test0,ddd_test2,tolerance = .000001)
  
  brts <- 1:10
  pars1 <- c(0.8,0,40)
  pars2 <- c(2,1,1,0,1)
  missnumspec <- 0
  result1 <- DDD::bd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec)
  result2 <- DDD::dd_loglik(pars1 = pars1,pars2 = c(min(1000,10 * (length(brts) + missnumspec)),1,pars2[2:5]),brts = brts,missnumspec = 0)
  testthat::expect_equal(result1,result2)
})           

context("test_DDD_KI_conditioning")

test_that("conditioning_DDD_KI works",
{          
  skip_if(Sys.getenv("CI") == "", message = "Run only on CI")
    ts <- seq(-9,-1,2);
    p1 <- rep(0,5);
    p2 <- rep(0,5);
    pars1_list <- list(c(0.5,0.4,Inf),c(0,0,Inf))
    reltol <- 1e-8
    abstol <- 1e-8

    lx_list <- list(200,200)
    for(i in 1:5)
    {
      brts_k_list <- list(rbind(c(-10,ts[i],0),c(2,1,1)),rbind(c(ts[i],0),c(1,1)))
      p1[i] <- DDD:::dd_multiple_KI_logliknorm(brts_k_list = brts_k_list,
                                               pars1_list = pars1_list,
                                               pars2 = c(200,1,5,NA,1,2,3),
                                               loglik = 0,
                                               lx_list = lx_list,
                                               reltol = reltol,
                                               abstol = abstol,
                                               methode = 'ode45')
      p2[i] <- DDD:::dd_KI_logliknorm(brts_k_list = brts_k_list,
                                      pars1_list = pars1_list,
                                      loglik = 0,
                                      cond = 5,
                                      ddep = 1,
                                      lx_list = lx_list,
                                      reltol = reltol,
                                      abstol = abstol,
                                      methode = 'ode45')
      #print(exp(p1[i]))
      #print(exp(p2[i]))
    }
    testthat::expect_equal(p1,p2,tolerance = 1e-4)


  pars1_list <- list(c(0.4,0.1,30),c(0.2,0.1,20))
  pars2 <- c(500,1,5,NA,1,2,3)
  lx_list <- list(pars2[1],pars2[1])
  brts_k_list <- list(rbind(sort(c(-10:-6,-3,-1,0)),c(2,3,4,5,6,5,6,6)),rbind(c(-3,-2,0),c(1,2,2)))
  logliknorm1 <- DDD:::dd_KI_logliknorm(brts_k_list = brts_k_list,
                                        pars1_list = pars1_list,
                                        loglik = 0,
                                        cond = 5,
                                        ddep = 1,
                                        lx_list = lx_list,
                                        reltol = reltol,
                                        abstol = abstol,
                                        methode = 'ode45')
  logliknorm2 <- DDD:::dd_multiple_KI_logliknorm(brts_k_list = brts_k_list,
                                                 pars1_list = pars1_list,
                                                 pars2 = pars2,
                                                 loglik = 0,
                                                 lx_list = lx_list,
                                                 reltol = reltol,
                                                 abstol = abstol,
                                                 methode = 'ode45')
  testthat::expect_equal(as.numeric(logliknorm1),logliknorm2,tolerance = 1e-4)
  
  pars1_list <- list(c(0.4,0.1,20),c(0.2,0.1,20))
  pars2 <- c(500,1,5,NA,1,2,3)
  lx_list <- list(pars2[1],pars2[1])
  brts_k_list <- list(rbind(c(-10:-6,-3,-1,0),c(2,3,4,5,6,5,6,6)),rbind(c(-3,-2,0),c(1,2,2)))
  logliknorm1 <- DDD:::dd_KI_logliknorm(brts_k_list = brts_k_list,
                                        pars1_list = pars1_list,
                                        loglik = 0,
                                        cond = 5,
                                        ddep = 1,
                                        lx_list = lx_list,
                                        reltol = reltol,
                                        abstol = abstol,
                                        methode = 'ode45')
  logliknorm2 <- DDD:::dd_multiple_KI_logliknorm(brts_k_list = brts_k_list,
                                                 pars1_list = pars1_list,
                                                 pars2 = pars2,
                                                 loglik = 0,
                                                 lx_list = lx_list,
                                                 reltol = reltol,
                                                 abstol = abstol,
                                                 methode = 'ode45')
  testthat::expect_equal(as.numeric(logliknorm1),logliknorm2,tolerance = 1e-4)
})