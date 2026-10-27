# Tests for bugs found in the review of develop (fae24c7) that the candidate patch does NOT
# fix. Each test states the intended behaviour and is skipped until the bug is fixed;
# remove the skip() call when it is. See ddd_test_report_develop_2026-10-26.md.

brts_9 <- c(10, 9, 7.5, 6, 4.2, 3, 2, 1.1, 0.5)

test_that("dd_loglik with cond = 3 and a large carrying capacity returns a number (report A.8)", {
  skip("Known bug: 'a' is computationally singular for roughly K * la / (la - mu) > 22-25")
  # The analytical cond = 3 normalisation uses a matrix inverse that becomes singular
  # for a larger state space; the error is not caught, so dd_ML(cond = 3) aborts.
  loglik <- DDD::dd_loglik(
    pars1 = c(0.8, 0.1, 25), pars2 = c(100, 1, 3, 1, 0, 2),
    brts = brts_9, missnumspec = 0
  )
  expect_true(is.numeric(loglik))
  expect_false(is.na(loglik))
})

test_that("rng_respecting_sample works with exactly one non-zero probability (report A.9)", {
  skip("Known bug: base::sample treats a single number as 1:x")
  set.seed(1)
  draws <- DDD::rng_respecting_sample(c(5, 6), size = 3, replace = TRUE, prob = c(1, 0))
  expect_equal(draws, c(5, 5, 5))
})

test_that("dd_KI_loglik with expm failure does not return a wrong finite value (report A.6)", {
  skip("Known bug: the try() fallback in dd_loglik_M is rep(0, length(p)) of length 1")
  # ddmodel 3, cond 4: the matrix exponential at res = 200 hits expm's mxrej limit. The
  # value at res = 100 (-100.0267) is the converged one; at res = 200 the fallback gives
  # -100.408.
  brtsM <- c(25.2, 24.6, 24.0, 22.5, 21.7, 20.4, 19.9, 19.7, 18.8, 17.1, 15.8, 11.8, 9.7, 8.9, 5.7, 5.2)
  brtsS <- c(9.6, 8.6, 7.4, 4.9, 2.5)
  pars1 <- c(0.25, 0.12, 25.51, 1.0, 0.16, 8.61, 9.8)
  loglik <- sapply(c(100, 200), function(res) {
    suppressWarnings(DDD::dd_KI_loglik(
      pars1 = pars1, pars2 = c(res, 3, 4, 18.8, 0, 2, 1),
      brtsM = brtsM, brtsS = brtsS, missnumspec = 0, methode = "analytical"
    ))
  })
  expect_equal(loglik[2], loglik[1], tolerance = 1e-4)
})

test_that("dd_sim works for ddmodels 1.4 and 1.5 (report A.3)", {
  skip("Known bug: dd_sim stops with \"object 'laN' not found\" for ddmodel 1.4 and 1.5")
  for (ddmodel in c(1.4, 1.5)) {
    set.seed(7)
    sim <- DDD::dd_sim(pars = c(0.8, 0.1, 25), age = 8, ddmodel = ddmodel)
    expect_true(length(sim$brts) >= 1, info = sprintf("ddmodel = %s", ddmodel))
  }
})

test_that("dd_loglik with cond = 3 and an ODE solver returns for ddmodels 2 and 3 (report A.2)", {
  skip("Known issue: after the logtrafo fix the ODE normalisation does not return within 90 s for ddmodels 2 and 3")
  for (ddmodel in c(2, 3)) {
    pars1 <- c(0.8, 0.1, 20)
    pars2 <- c(100, ddmodel, 3, 1, 0, 2)
    analytical <- DDD::dd_loglik(pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0)
    numerical <- DDD::dd_loglik(
      pars1 = pars1, pars2 = pars2, brts = brts_9, missnumspec = 0,
      methode = "odeint::runge_kutta_cash_karp54"
    )
    expect_equal(as.numeric(numerical), as.numeric(analytical), tolerance = 1e-4,
      info = sprintf("ddmodel = %s", ddmodel))
  }
})

test_that("dd_SR_loglik with cond = 3 works with an ODE solver (report B.10)", {
  skip("Known bug: the cond = 3 branch of dd_SR_loglik1 calls deSolve::ode directly")
  pars <- c(0.8, 0.1, 20)
  pars2 <- c(100, 1, 3, 1, 0, 2)
  analytical <- DDD::dd_SR_loglik(pars1 = c(pars, pars, 5), pars2 = pars2, brts = brts_9, missnumspec = 0)
  numerical <- DDD::dd_SR_loglik(
    pars1 = c(pars, pars, 5), pars2 = pars2, brts = brts_9, missnumspec = 0,
    methode = "odeint::runge_kutta_cash_karp54"
  )
  expect_equal(as.numeric(numerical), as.numeric(analytical), tolerance = 1e-4)
})
