library(DDD)


methods = c(
  'analytical',
  'ode45',
  'lsoda',
  'odeint::runge_kutta_cash_karp54',
  'odeint::runge_kutta_fehlberg78',    # default
  'odeint::runge_kutta_dopri5',
  'odeint::bulirsch_stoer'
)


test_ode = function(rep = 1) {
  # parameter set from test suite:
  pars1 = c(0.8,0.1,40) 
  pars2 = c(100,1,3,0.1,0.1,2)  # pars2[3] = 1: fw; = 3: bw;
  brts = 1:30
  missnumspec = 0

  for (m in methods) 
  {
#    print(system.time(
      for (i in 1:rep) {
        r <- DDD:::dd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,methode = m)
      }
#    ))
    print(paste(m, r), digits=17)
    print('-----------')
  }
}


test_KI = function(imeth) {
  pars <- c(0.4, 0.1, 0.6, 0.2)
  brts <- list(10:6, 3:2)
  high_k <- 1e7
  pars1 <- c(pars[1], pars[2], high_k, pars[3], pars[4], high_k, brts[[2]][1])
  pars2 <- c(200,1,0,brts[[2]][1],0,2,1.5)
  brtsM <- brts[[1]]
  brtsS <- brts[[2]][-1]
  
  pars2[3] <- 1
  print(methods[imeth])
  ddd_test <- DDD::dd_KI_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brtsM = brtsM,
    brtsS = brtsS,
    missnumspec = 0,
    methode = methods[imeth]
  ); 
  
  print(ddd_test)
}
