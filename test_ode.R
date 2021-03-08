library(DDD)


test_ode = function() {
  # parameter set from test suite:
  pars1 = c(0.8,0.1,40) 
  pars2 = c(100,1,1,0,0,2)
  brts = 1:30
  missnumspec = 0

  methods = c('lsoda',
              'ode45',
              'odeint::runge_kutta_cash_karp54',
              'odeint::runge_kutta_fehlberg78',    # default
              'odeint::runge_kutta_dopri5',
              'odeint::bulirsch_stoer')

  for (m in methods) 
  {
    r <- DDD:::dd_loglik(pars1 = pars1,pars2 = pars2,brts = brts,missnumspec = missnumspec,methode = m)
    print(paste(m, r), digits=17)
  }
}
