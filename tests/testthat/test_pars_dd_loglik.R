# Case tree
set.seed(26071948) # jjsepkoski
phylo <- dd_sim(
  pars = c(0.8, 0.1, 20, 1),
  age = 20,
  ddmodel = 5
)$tes
brts <- ape::branching.times(phylo)

# Shortcut function
run_dd_loglik <- function(ddmodel, lambda_0 = 0.8, mu_0 = 0.1, K = 20, r = 1, verbose = FALSE) {
  if (both_rates_vary(ddmodel)) {
    pars1 <- c(lambda_0, mu_0, K, r)
  } else {
    pars1 <- c(lambda_0, mu_0, K)
  }
  pars2 <- c("lx" = 50, "ddmodel" = ddmodel, "cond" = 1, "btorph" = 0, "verbose" = verbose, "soc" = 2)
  dd_loglik(
    pars1 = pars1,
    pars2 = pars2,
    brts = brts, # defined above
    missnumspec = 0,
    methode = "odeint::runge_kutta_cash_karp54"
  )
}

test_that("all DD models return a likelihood", {
  expect_silent(logL_dd1 <- run_dd_loglik(ddmodel = 1))
  expect_true(logL_dd1 < 0 && is.finite(logL_dd1))
  expect_silent(logL_dd2 <- run_dd_loglik(ddmodel = 2))
  expect_true(logL_dd2 < 0 && is.finite(logL_dd2))
  expect_silent(logL_dd3 <- run_dd_loglik(ddmodel = 3))
  expect_true(logL_dd3 < 0 && is.finite(logL_dd3))
  expect_silent(logL_dd4 <- run_dd_loglik(ddmodel = 4))
  expect_true(logL_dd4 < 0 && is.finite(logL_dd4))
  expect_silent(logL_dd5 <- run_dd_loglik(ddmodel = 5))
  expect_true(logL_dd5 < 0 && is.finite(logL_dd5))
  expect_silent(logL_dd6 <- run_dd_loglik(ddmodel = 6))
  expect_true(logL_dd6 < 0 && is.finite(logL_dd6))
  expect_silent(logL_dd7 <- run_dd_loglik(ddmodel = 7))
  expect_true(logL_dd7 < 0 && is.finite(logL_dd7))
  expect_silent(logL_dd8 <- run_dd_loglik(ddmodel = 8))
  expect_true(logL_dd8 < 0 && is.finite(logL_dd8))
  expect_silent(logL_dd9 <- run_dd_loglik(ddmodel = 9))
  expect_true(logL_dd9 < 0 && is.finite(logL_dd9))
  expect_silent(logL_dd10 <- run_dd_loglik(ddmodel = 10))
  expect_true(logL_dd10 < 0 && is.finite(logL_dd10))
  expect_silent(logL_dd11 <- run_dd_loglik(ddmodel = 11))
  expect_true(logL_dd11 < 0 && is.finite(logL_dd11))
  expect_silent(logL_dd12 <- run_dd_loglik(ddmodel = 12))
  expect_true(logL_dd12 < 0 && is.finite(logL_dd12))
  expect_silent(logL_dd13 <- run_dd_loglik(ddmodel = 13))
  expect_true(logL_dd13 < 0 && is.finite(logL_dd13))
  expect_silent(logL_dd14 <- run_dd_loglik(ddmodel = 14))
  expect_true(logL_dd14 < 0 && is.finite(logL_dd14))
  expect_silent(logL_dd15 <- run_dd_loglik(ddmodel = 15))
  expect_true(logL_dd15 < 0 && is.finite(logL_dd15))
})

test_that("forbidden parameter values are handled properly", {
  expect_equal(run_dd_loglik(ddmodel = 1, lambda_0 = 10000), -Inf)
  # Negative parameters are not allowed
  expect_equal(run_dd_loglik(ddmodel = 1, lambda_0 = -1), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 5, mu_0 = -1), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 1, K = -1), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 15, r = -1), -Inf)
  # lambda0 = mu0 not allowed
  expect_equal(run_dd_loglik(ddmodel = 1, lambda_0 = 0.8, mu_0 = 0.8), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 1, lambda_0 = 0, mu_0 = 0, K = 0), -Inf)
  # Linear-DD speciation must have K' > N
  expect_equal(run_dd_loglik(ddmodel = 1, K = 14), -Inf) 
  expect_equal(run_dd_loglik(ddmodel = 5, K = 7), -Inf)
  # Other values of K < N but respecting K' > N should return a result
  expect_silent(run_dd_loglik(ddmodel = 5, K = 8))
  # Exponential-DD extinction must have mu0 > 0 or result is NaN
  expect_equal(run_dd_loglik(ddmodel = 4, mu_0 = 0), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 6, mu_0 = 0), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 7, mu_0 = 0), -Inf)
  # Exponential-DD speciation must have either mu0 > 0 or r > 0 or result is NaN
  expect_equal(run_dd_loglik(ddmodel = 2, mu_0 = 0), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 8, mu_0 = 0, r = 0), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 14, mu_0 = 0, r = 0), -Inf)
  # Exponential(alternative)-DD extinction must have mu0 > 0 or result is NaN
  expect_equal(run_dd_loglik(ddmodel = 10, mu_0 = 0), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 11, mu_0 = 0), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 12, mu_0 = 0), -Inf)
  expect_equal(run_dd_loglik(ddmodel = 15, mu_0 = 0, r = 0), -Inf)
})