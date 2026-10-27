# dd_loglik: the analytical method (matrix exponential) and the ODE solver must agree.
#
# On develop (fae24c7):
# - ddmodels 2.1, 2.2 and 2.3 with an ODE solver and stem age (soc = 1) stop with
#   "missing value where TRUE/FALSE needed": in dd_loglik_rhs_precomp the speciation
#   rate at n = 0 is la * 0^(-x) = Inf, which is multiplied by zero padding and gives NaN.
# - cond = 3 with an ODE solver stops with "object 'logtrafo' not found": logtrafo and
#   lx_old are only created inside if(cond != 3) in dd_loglik1.

brts_9 <- c(10, 9, 7.5, 6, 4.2, 3, 2, 1.1, 0.5)
ode <- "odeint::runge_kutta_cash_karp54"

pars_for <- function(ddmodel) {
  if (ddmodel == 2.3) return(c(0.8, 0.1, 0.5))
  if (ddmodel == 5) return(c(0.8, 0.1, 20, 1))
  c(0.8, 0.1, 20)
}

expect_ode_equals_analytical <- function(ddmodel, cond, soc, tolerance = 1e-5) {
  pars1 <- pars_for(ddmodel)
  pars2 <- c(60, ddmodel, cond, 1, 0, soc)
  analytical <- DDD::dd_loglik(
    pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0
  )
  numerical <- DDD::dd_loglik(
    pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0, methode = ode
  )
  expect_true(is.finite(analytical), info = sprintf("ddmodel = %s, cond = %d, soc = %d", ddmodel, cond, soc))
  expect_equal(
    as.numeric(numerical), as.numeric(analytical), tolerance = tolerance,
    info = sprintf("ddmodel = %s, cond = %d, soc = %d", ddmodel, cond, soc)
  )
}

test_that("dd_loglik: ODE equals analytical for crown and stem age, ddmodels 1, 1.3, 2, 3, 4, 4.1, 4.2, 5", {
  skip_on_cran()
  for (ddmodel in c(1, 1.3, 2, 3, 4, 4.1, 4.2, 5)) {
    for (cond in 0:1) {
      for (soc in 1:2) {
        expect_ode_equals_analytical(ddmodel, cond, soc)
      }
    }
  }
})

test_that("dd_loglik: ODE equals analytical for ddmodels 2.1, 2.2, 2.3 with crown age", {
  skip_on_cran()
  for (ddmodel in c(2.1, 2.2, 2.3)) {
    for (cond in 0:1) {
      expect_ode_equals_analytical(ddmodel, cond, soc = 2)
    }
  }
})

test_that("dd_loglik: ODE equals analytical for ddmodels 2.1, 2.2, 2.3 with stem age", {
  skip_on_cran()
  for (ddmodel in c(2.1, 2.2, 2.3)) {
    for (cond in 0:1) {
      expect_ode_equals_analytical(ddmodel, cond, soc = 1)
    }
  }
})

test_that("dd_loglik: ODE equals analytical for cond = 3 (ddmodel 1)", {
  skip_on_cran()
  pars1 <- c(0.8, 0.1, 20)
  pars2 <- c(100, 1, 3, 1, 0, 2)
  analytical <- DDD::dd_loglik(
    pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0
  )
  numerical <- DDD::dd_loglik(
    pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0, methode = ode
  )
  # the normalisation by the ODE is integrated to a large finite time, so the two
  # methods differ by about 1e-4 (3e-6 relative)
  expect_equal(as.numeric(numerical), as.numeric(analytical), tolerance = 1e-4)
})

test_that("dd_loglik: ODE equals analytical with missing species", {
  skip_on_cran()
  pars1 <- c(0.8, 0.1, 20)
  for (cond in 0:2) {
    pars2 <- c(100, 1, cond, 1, 0, 2)
    analytical <- DDD::dd_loglik(
      pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 2
    )
    numerical <- DDD::dd_loglik(
      pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 2, methode = ode
    )
    expect_equal(as.numeric(numerical), as.numeric(analytical), tolerance = 1e-5,
      info = sprintf("cond = %d", cond))
  }
})
