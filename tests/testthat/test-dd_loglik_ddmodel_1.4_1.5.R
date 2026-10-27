# ddmodels 1.4 and 1.5 are documented in dd_loglik and dd_ML. Their analytical and ODE
# likelihoods must agree.
#
# On develop (fae24c7):
# - lambdamu() has no branch for ddmodel 1.4 or 1.5, so the analytical likelihood is the
#   constant-rate likelihood and does not depend on K (-41.07678 for every K).
# - dd_loglik_rhs has no muvec for ddmodel 1.4: the ODE solver stops with
#   "object 'muvec' not found".
#
# For ddmodel 1.4 with cond = 1 the two methods still differ by about 0.006 at res = 100
# (the speciation rate la * n / (n + K) decays slowly, so the conditioning sum over the
# truncated state space is not converged); cond = 0 and cond = 2 agree.

brts_9 <- c(10, 9, 7.5, 6, 4.2, 3, 2, 1.1, 0.5)
ode <- "odeint::runge_kutta_cash_karp54"

test_that("ddmodel 1.4: analytical equals ODE for different carrying capacities", {
  skip_on_cran()
  for (K in c(5, 20, 100)) {
    for (cond in c(0, 2)) {
      pars2 <- c(100, 1.4, cond, 1, 0, 2)
      analytical <- DDD::dd_loglik(
        pars1 = c(0.8, 0.1, K), pars2 = pars2, brts = brts_9, missnumspec = 0
      )
      numerical <- DDD::dd_loglik(
        pars1 = c(0.8, 0.1, K), pars2 = pars2, brts = brts_9, missnumspec = 0, methode = ode
      )
      expect_equal(as.numeric(numerical), as.numeric(analytical), tolerance = 1e-5,
        info = sprintf("K = %g, cond = %d", K, cond))
    }
  }
})

test_that("ddmodel 1.4: the analytical likelihood depends on K", {
  skip_on_cran()
  loglik <- sapply(c(5, 20, 100), function(K) {
    DDD::dd_loglik(
      pars1 = c(0.8, 0.1, K), pars2 = c(100, 1.4, 1, 1, 0, 2),
      brts = brts_9, missnumspec = 0
    )
  })
  expect_true(all(is.finite(loglik)))
  expect_true(min(abs(diff(sort(loglik)))) > 0.1)
})

test_that("ddmodel 1.5: analytical equals ODE", {
  skip_on_cran()
  for (K in c(25, 40)) {
    for (cond in 0:2) {
      pars2 <- c(100, 1.5, cond, 1, 0, 2)
      analytical <- DDD::dd_loglik(
        pars1 = c(0.8, 0.1, K), pars2 = pars2, brts = brts_9, missnumspec = 0
      )
      numerical <- DDD::dd_loglik(
        pars1 = c(0.8, 0.1, K), pars2 = pars2, brts = brts_9, missnumspec = 0, methode = ode
      )
      expect_equal(as.numeric(numerical), as.numeric(analytical), tolerance = 1e-5,
        info = sprintf("K = %g, cond = %d", K, cond))
    }
  }
})
