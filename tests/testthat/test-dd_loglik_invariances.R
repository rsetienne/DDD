# Invariances and limiting cases of dd_loglik and bd_loglik. These pass on develop
# (fae24c7); they guard against regressions in the fixes tested in the other files.

brts_9 <- c(10, 9, 7.5, 6, 4.2, 3, 2, 1.1, 0.5)

test_that("dd_loglik does not depend on the order or the sign of brts", {
  pars1 <- c(0.8, 0.1, 20)
  pars2 <- c(100, 1, 1, 0, 0, 2)
  reference <- DDD::dd_loglik(pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0)
  expect_equal(
    DDD::dd_loglik(pars1 = pars1, pars2 = pars2, brts = rev(brts_9), missnumspec = 0),
    reference
  )
  expect_equal(
    DDD::dd_loglik(pars1 = pars1, pars2 = pars2, brts = -brts_9, missnumspec = 0),
    reference
  )
})

test_that("dd_loglik: btorph = 0 versus 1 differs by lgamma of the number of tips", {
  pars1 <- c(0.8, 0.1, 20)
  for (soc in 1:2) {
    n_tips <- length(brts_9) + soc - 1
    ll0 <- DDD::dd_loglik(pars1 = pars1, pars2 = c(100, 1, 1, 0, 0, soc), brts = brts_9, missnumspec = 0)
    ll1 <- DDD::dd_loglik(pars1 = pars1, pars2 = c(100, 1, 1, 1, 0, soc), brts = brts_9, missnumspec = 0)
    expect_equal(as.numeric(ll0 - ll1), lgamma(n_tips), info = sprintf("soc = %d", soc))
  }
})

test_that("dd_loglik with an infinite carrying capacity equals bd_loglik", {
  for (soc in 1:2) {
    for (cond in 0:2) {
      for (missnumspec in c(0, 2)) {
        expected <- DDD::bd_loglik(
          pars1 = c(0.5, 0.1), pars2 = c(0, cond, 1, 0, soc),
          brts = brts_9, missnumspec = missnumspec
        )
        observed <- DDD::dd_loglik(
          pars1 = c(0.5, 0.1, Inf), pars2 = c(100, 1, cond, 1, 0, soc),
          brts = brts_9, missnumspec = missnumspec
        )
        expect_equal(as.numeric(observed), as.numeric(expected),
          info = sprintf("soc = %d, cond = %d, missnumspec = %d", soc, cond, missnumspec))
      }
    }
  }
})

test_that("bd_loglik does not depend on the order or the sign of brts", {
  pars1 <- c(0.5, 0.1)
  pars2 <- c(0, 1, 0, 0, 2)
  reference <- DDD::bd_loglik(pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0)
  expect_equal(DDD::bd_loglik(pars1 = pars1, pars2 = pars2, brts = rev(brts_9), missnumspec = 0), reference)
  expect_equal(DDD::bd_loglik(pars1 = pars1, pars2 = pars2, brts = -brts_9, missnumspec = 0), reference)
})

test_that("dd_SR_loglik with cond = 3 and crown age equals dd_loglik without a shift", {
  skip_on_cran()
  pars <- c(0.8, 0.1, 20)
  pars2 <- c(100, 1, 3, 1, 0, 2)
  expected <- DDD::dd_loglik(pars1 = pars, pars2 = pars2, brts = brts_9, missnumspec = 0)
  observed <- DDD::dd_SR_loglik(pars1 = c(pars, pars, 5), pars2 = pars2, brts = brts_9, missnumspec = 0)
  expect_equal(as.numeric(observed), as.numeric(expected), tolerance = 1e-5)
})
