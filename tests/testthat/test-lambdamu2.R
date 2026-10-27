# lambdamu2 builds the speciation and extinction rate matrices used for the two-group
# (KI and MS) conditioning. For the models it supports it must give the same rates as
# lambdamu.
#
# On develop (fae24c7) lambdamu2 has muvec = (nxt + n0)^x for ddmodels 4, 4.1 and 4.2,
# which lacks the factor mu (lambdamu has muM * (n + n0)^x): the extinction rates are
# about 1 / mu = 10 times too large for mu = 0.1. As a result dd_KI_loglik with
# ddmodel = 4 and cond = 4 gives -Inf with the analytical method.

test_that("lambdamu2 gives the same rates as lambdamu", {
  n <- 0:3
  nxt <- outer(n, n, "+")
  pars <- c(0.8, 0.1, 20)
  for (ddmodel in c(1, 1.3, 2, 3, 4, 4.1, 4.2)) {
    rates <- DDD:::lambdamu(as.numeric(nxt), c(pars, 0), ddmodel)
    rates2 <- DDD:::lambdamu2(nxt, pars, ddmodel)
    expect_equal(as.numeric(rates2[[1]]), as.numeric(rates[[1]]),
      info = sprintf("speciation rates, ddmodel = %s", ddmodel))
    expect_equal(as.numeric(rates2[[2]]), as.numeric(rates[[2]]),
      info = sprintf("extinction rates, ddmodel = %s", ddmodel))
  }
})

test_that("dd_KI_loglik with ddmodel 4 and cond = 4: analytical equals ODE", {
  skip_on_cran()
  brtsM <- c(25.2, 24.6, 24.0, 22.5, 21.7, 20.4, 19.9, 19.7, 18.8, 17.1, 15.8, 11.8, 9.7, 8.9, 5.7, 5.2)
  brtsS <- c(9.6, 8.6, 7.4, 4.9, 2.5)
  pars1 <- c(0.25, 0.12, 25.51, 1.0, 0.16, 8.61, 9.8)
  pars2 <- c(200, 4, 4, 18.8, 0, 2, 1)
  analytical <- DDD::dd_KI_loglik(
    pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS,
    missnumspec = 0, methode = "analytical"
  )
  numerical <- DDD::dd_KI_loglik(
    pars1 = pars1, pars2 = pars2, brtsM = brtsM, brtsS = brtsS,
    missnumspec = 0, methode = "odeint::runge_kutta_cash_karp54"
  )
  expect_true(is.finite(analytical))
  expect_equal(as.numeric(analytical), as.numeric(numerical), tolerance = 1e-5)
})
