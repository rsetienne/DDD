# dd_SR_loglik with identical parameters before and after the shift must give the
# same likelihood as dd_loglik, for crown age (soc = 2) and for stem age (soc = 1).
#
# On develop (fae24c7) the soc = 1 comparisons fail: in dd_SR_loglik1 and dd_SR_loglik2
# the term (cond != 3 & soc == 2) * log(probs[...]) drops the log-probability for stem
# age (dd_loglik uses (cond != 3 | soc == 1)), and probsc[S + missnumspec - 1] should be
# probsc[S + missnumspec - soc + 1] for cond = 2.

brts_9 <- c(10, 9, 7.5, 6, 4.2, 3, 2, 1.1, 0.5)

test_that("dd_SR_loglik without a shift in parameters equals dd_loglik (analytical)", {
  skip_on_cran()
  pars <- c(0.8, 0.1, 20)
  for (ddmodel in c(1, 1.3, 2, 3, 4)) {
    for (cond in 0:2) {
      for (soc in 1:2) {
        pars2 <- c(100, ddmodel, cond, 1, 0, soc)
        expected <- DDD::dd_loglik(
          pars1 = pars, pars2 = pars2, brts = brts_9, missnumspec = 0
        )
        observed <- DDD::dd_SR_loglik(
          pars1 = c(pars, pars, 5), pars2 = pars2, brts = brts_9, missnumspec = 0
        )
        expect_equal(
          as.numeric(observed), as.numeric(expected), tolerance = 1e-5,
          info = sprintf("ddmodel = %s, cond = %d, soc = %d", ddmodel, cond, soc)
        )
      }
    }
  }
})

test_that("dd_SR_loglik without a shift in parameters equals dd_loglik (missing species)", {
  skip_on_cran()
  pars <- c(0.8, 0.1, 20)
  for (cond in 0:2) {
    for (soc in 1:2) {
      pars2 <- c(100, 1, cond, 1, 0, soc)
      expected <- DDD::dd_loglik(
        pars1 = pars, pars2 = pars2, brts = brts_9, missnumspec = 2
      )
      observed <- DDD::dd_SR_loglik(
        pars1 = c(pars, pars, 5), pars2 = pars2, brts = brts_9, missnumspec = 2
      )
      expect_equal(
        as.numeric(observed), as.numeric(expected), tolerance = 1e-5,
        info = sprintf("cond = %d, soc = %d, missnumspec = 2", cond, soc)
      )
    }
  }
})

test_that("dd_SR_loglik with an ODE solver equals dd_loglik for stem age", {
  skip_on_cran()
  pars <- c(0.8, 0.1, 20)
  for (cond in 0:2) {
    pars2 <- c(100, 1, cond, 1, 0, 1)
    expected <- DDD::dd_loglik(
      pars1 = pars, pars2 = pars2, brts = brts_9, missnumspec = 0
    )
    observed <- DDD::dd_SR_loglik(
      pars1 = c(pars, pars, 5), pars2 = pars2, brts = brts_9, missnumspec = 0,
      methode = "odeint::runge_kutta_cash_karp54"
    )
    expect_equal(
      as.numeric(observed), as.numeric(expected), tolerance = 1e-4,
      info = sprintf("cond = %d, soc = 1", cond)
    )
  }
})

test_that("dd_SR_loglik analytical and ODE agree when the parameters shift", {
  skip_on_cran()
  pars1 <- c(0.8, 0.1, 20, 0.4, 0.2, 12, 5)
  for (cond in 0:2) {
    for (soc in 1:2) {
      pars2 <- c(100, 1, cond, 1, 0, soc)
      analytical <- DDD::dd_SR_loglik(
        pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0
      )
      ode <- DDD::dd_SR_loglik(
        pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0,
        methode = "odeint::runge_kutta_cash_karp54"
      )
      expect_equal(
        as.numeric(analytical), as.numeric(ode), tolerance = 1e-5,
        info = sprintf("cond = %d, soc = %d", cond, soc)
      )
    }
  }
})
