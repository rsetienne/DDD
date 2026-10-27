# With cond = 3 the likelihood is conditioned on the total number of extant taxa only and is
# computed by integrating over all possible stem ages. The first branching time is therefore
# always the crown age: pars2[6] (soc) must have no effect when cond = 3.
#
# On develop (fae24c7) dd_loglik and bd_loglik do this, but dd_SR_loglik does not: with
# soc = 1 it reads the first branching time as a stem age (-41.37 instead of -32.30 for
# the example below).

brts_9 <- c(10, 9, 7.5, 6, 4.2, 3, 2, 1.1, 0.5)

test_that("dd_loglik: soc has no effect when cond = 3", {
  pars1 <- c(0.8, 0.1, 20)
  for (ddmodel in c(1, 2, 3, 4)) {
    crown <- DDD::dd_loglik(
      pars1 = pars1, pars2 = c(100, ddmodel, 3, 1, 0, 2), brts = brts_9, missnumspec = 0
    )
    stem <- DDD::dd_loglik(
      pars1 = pars1, pars2 = c(100, ddmodel, 3, 1, 0, 1), brts = brts_9, missnumspec = 0
    )
    expect_true(is.finite(crown), info = sprintf("ddmodel = %s", ddmodel))
    expect_equal(as.numeric(stem), as.numeric(crown), info = sprintf("ddmodel = %s", ddmodel))
  }
})

test_that("dd_loglik with an ODE solver: soc has no effect when cond = 3", {
  skip_on_cran()
  pars1 <- c(0.8, 0.1, 20)
  ode <- "odeint::runge_kutta_cash_karp54"
  crown <- DDD::dd_loglik(
    pars1 = pars1, pars2 = c(100, 1, 3, 1, 0, 2), brts = brts_9, missnumspec = 0, methode = ode
  )
  stem <- DDD::dd_loglik(
    pars1 = pars1, pars2 = c(100, 1, 3, 1, 0, 1), brts = brts_9, missnumspec = 0, methode = ode
  )
  expect_equal(as.numeric(stem), as.numeric(crown))
})

test_that("dd_loglik: soc has no effect when cond = 3 with missing species", {
  pars1 <- c(0.8, 0.1, 20)
  crown <- DDD::dd_loglik(
    pars1 = pars1, pars2 = c(100, 1, 3, 1, 0, 2), brts = brts_9, missnumspec = 2
  )
  stem <- DDD::dd_loglik(
    pars1 = pars1, pars2 = c(100, 1, 3, 1, 0, 1), brts = brts_9, missnumspec = 2
  )
  expect_equal(as.numeric(stem), as.numeric(crown))
})

test_that("dd_loglik: soc does have an effect for the other conditionings", {
  pars1 <- c(0.8, 0.1, 20)
  for (cond in 0:2) {
    crown <- DDD::dd_loglik(
      pars1 = pars1, pars2 = c(100, 1, cond, 1, 0, 2), brts = brts_9, missnumspec = 0
    )
    stem <- DDD::dd_loglik(
      pars1 = pars1, pars2 = c(100, 1, cond, 1, 0, 1), brts = brts_9, missnumspec = 0
    )
    expect_gt(abs(stem - crown), 0.1, label = sprintf("cond = %d", cond))
  }
})

test_that("bd_loglik: soc has no effect when cond = 3", {
  crown <- DDD::bd_loglik(
    pars1 = c(0.5, 0.1), pars2 = c(0, 3, 1, 0, 2), brts = brts_9, missnumspec = 0
  )
  stem <- DDD::bd_loglik(
    pars1 = c(0.5, 0.1), pars2 = c(0, 3, 1, 0, 1), brts = brts_9, missnumspec = 0
  )
  expect_equal(as.numeric(stem), as.numeric(crown))
})

test_that("dd_SR_loglik: soc has no effect when cond = 3 and equals dd_loglik without a shift", {
  skip_on_cran()
  pars <- c(0.8, 0.1, 20)
  for (ddmodel in c(1, 2, 3)) {
    reference <- DDD::dd_loglik(
      pars1 = pars, pars2 = c(100, ddmodel, 3, 1, 0, 2), brts = brts_9, missnumspec = 0
    )
    for (soc in 1:2) {
      observed <- DDD::dd_SR_loglik(
        pars1 = c(pars, pars, 5), pars2 = c(100, ddmodel, 3, 1, 0, soc),
        brts = brts_9, missnumspec = 0
      )
      expect_equal(as.numeric(observed), as.numeric(reference), tolerance = 1e-5,
        info = sprintf("ddmodel = %s, soc = %d", ddmodel, soc))
    }
  }
})

test_that("dd_SR_loglik: soc has no effect when cond = 3 and the parameters shift", {
  skip_on_cran()
  pars1 <- c(0.8, 0.1, 20, 0.4, 0.2, 12, 5)
  crown <- DDD::dd_SR_loglik(
    pars1 = pars1, pars2 = c(100, 1, 3, 1, 0, 2), brts = brts_9, missnumspec = 0
  )
  stem <- DDD::dd_SR_loglik(
    pars1 = pars1, pars2 = c(100, 1, 3, 1, 0, 1), brts = brts_9, missnumspec = 0
  )
  expect_true(is.finite(crown))
  expect_equal(as.numeric(stem), as.numeric(crown))
})
