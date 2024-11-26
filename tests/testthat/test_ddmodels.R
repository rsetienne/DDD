context("test_ddmodels")

# Case 1.: 0 < r < Inf (or 0 < phi < 1)
pars_set1 <- c(
  "lambda_0" = 0.8,
  "mu_0" = 0.2,
  "K" = 20,
  "r" = 1/3 # corresponds to phi = 1/4
)
# Rates obtained on paper
exptd_rates_set1 <- list(
  "lambda_cst"    = function(N) 0.8, # not with phi != 1
  "mu_cst"        = function(N) 0.2, # not with phi != 0
  "lambda_lin"     = function(N) pmax(0.8 - 0.0225 * N, 0),
  "mu_lin"         = function(N) 0.2 + 0.0075 * N,
  "lambda_pow"     = function(N) pmax(0.8 * N ^ (-log(0.8 / 0.35) / log(20)), 0),
  "mu_pow"         = function(N) 0.2 * N ^ (log(1.75) / log(20)),
  "lambda_exp" = function(N) pmax(0.8 * (7 / 16) ^ (N / 20), 0),
  'mu_exp'     = function(N) 0.2 * (7 / 4) ^ (N / 20)
)

# Case 2.: r = 0 (phi = 0)
pars_set2 <- c(
  "lambda_0" = 0.8,
  "mu_0" = 0.2,
  "K" = 20,
  "r" = 0
)
exptd_rates_set2 <- list(
  "lambda_cst"    = function(N) rep(0.8, length(N)), # not with phi != 1
  "mu_cst"        = function(N) rep(0.2, length(N)), # not with phi != 0
  "lambda_lin"     = function(N) pmax(0.8 - 0.03 * N, 0),
  "mu_lin"         = function(N) rep(0.2, length(N)),
  "lambda_pow"     = function(N) pmax(0.8 * N ^ (-log(4) / log(20)), 0),
  "mu_pow"         = function(N) rep(0.2, length(N)),
  "lambda_exp" = function(N) pmax(0.8 * (1 / 4) ^ (N / 20), 0),
  'mu_exp'     = function(N) rep(0.2, length(N))
)
# Case 3.: r = Inf (phi = 1)
pars_set3 <- c(
  "lambda_0" = 0.8,
  "mu_0" = 0.2,
  "K" = 20,
  "r" = Inf
)
exptd_rates_set3 <- list(
  "lambda_cst"    = function(N) rep(0.8, length(N)), # not with phi != 1
  "mu_cst"        = function(N) rep(0.2, length(N)), # not with phi != 0
  "lambda_lin"     = function(N) rep(0.8, length(N)),
  "mu_lin"         = function(N) 0.2 + 0.03 * N,
  "lambda_pow"     = function(N) rep(0.8, length(N)),
  "mu_pow"         = function(N) 0.2 * N ^ (log(4) / log(20)),
  "lambda_exp" = function(N) rep(0.8, length(N)),
  'mu_exp'     = function(N) 0.2 * 4 ^ (N / 20)
)

# Declare test functions
## Match a ddmodel with the corresponding pair of speciation and extinction functions
match_exptd_rates <- function(ddmodel, exptd_rates_set, n_seq) {
  rates_ls <- switch(
    as.character(ddmodel),
    "1" = list("la_N" = exptd_rates_set$lambda_lin(n_seq), "mu_N" = exptd_rates_set$mu_cst(n_seq)),
    "2" = list("la_N" = exptd_rates_set$lambda_pow(n_seq), "mu_N" = exptd_rates_set$mu_cst(n_seq)),
    "3" = list("la_N" = exptd_rates_set$lambda_cst(n_seq), "mu_N" = exptd_rates_set$mu_lin(n_seq)),
    "4" = list("la_N" = exptd_rates_set$lambda_cst(n_seq), "mu_N" = exptd_rates_set$mu_pow(n_seq)),
    "5" = list("la_N" = exptd_rates_set$lambda_lin(n_seq), "mu_N" = exptd_rates_set$mu_lin(n_seq)),
    "6" = list("la_N" = exptd_rates_set$lambda_lin(n_seq), "mu_N" = exptd_rates_set$mu_pow(n_seq)),
    "7" = list("la_N" = exptd_rates_set$lambda_pow(n_seq), "mu_N" = exptd_rates_set$mu_pow(n_seq)),
    "8" = list("la_N" = exptd_rates_set$lambda_pow(n_seq), "mu_N" = exptd_rates_set$mu_lin(n_seq)),
    "9" = list("la_N" = exptd_rates_set$lambda_exp(n_seq), "mu_N" = exptd_rates_set$mu_cst(n_seq)),
    "10" = list("la_N" = exptd_rates_set$lambda_cst(n_seq), "mu_N" = exptd_rates_set$mu_exp(n_seq)),
    "11" = list("la_N" = exptd_rates_set$lambda_lin(n_seq), "mu_N" = exptd_rates_set$mu_exp(n_seq)),
    "12" = list("la_N" = exptd_rates_set$lambda_exp(n_seq), "mu_N" = exptd_rates_set$mu_exp(n_seq)),
    "13" = list("la_N" = exptd_rates_set$lambda_exp(n_seq), "mu_N" = exptd_rates_set$mu_lin(n_seq)),
    "14" = list("la_N" = exptd_rates_set$lambda_exp(n_seq), "mu_N" = exptd_rates_set$mu_pow(n_seq)),
    "15" = list("la_N" = exptd_rates_set$lambda_pow(n_seq), "mu_N" = exptd_rates_set$mu_exp(n_seq))
  )
  return(rates_ls)
}

## Test function; compare DDD output with rates obtained on paper
test_dd_loglik_rhs_precomp <- function(ddmodel, pars_set, exptd_rates_set) {
  # global variables
  N <- pars_set["K"]
  x <- rep(NA, 10) # length of the Q_n vector
  n_seq <- c(0, 0:(length(x) + 2 * N)) # based on internal code, not sure why
  lnn <- length(n_seq)
  
  exptd_rates <- match_exptd_rates(
    ddmodel = ddmodel,
    exptd_rates_set = exptd_rates_set,
    n_seq = n_seq
  )
  ddd_output <- dd_loglik_rhs_precomp(
    pars = c("pars" = pars_set, "k" = N, "ddmodel" = ddmodel), x = x
  )
  names(ddd_output) <- NULL
  ddd_rates <- list("la_N" = ddd_output[1:lnn], "mu_N" = ddd_output[(lnn + 1):(2 * lnn)])
  # Test
  cat(paste("Testing ddmodel =", ddmodel, "\n"))
  expect_equal(ddd_rates, exptd_rates)
}

test_lambdamu <- function(ddmodel, pars_set, exptd_rates_set) {
  # global variables
  N <- pars_set["K"]
  x <- rep(NA, 10) # length of the Q_n vector
  n_seq <- c(0, 0:(length(x) + 2 * N)) # based on internal code, not sure why
  lnn <- length(n_seq)
  
  exptd_rates <- match_exptd_rates(
    ddmodel = ddmodel,
    exptd_rates_set = exptd_rates_set,
    n_seq = n_seq
  )
  ddd_rates <- lambdamu(n = n_seq, pars = pars_set, ddep = ddmodel)
  names(ddd_rates) <- c("la_N", "mu_N")
  names(ddd_rates$la_N) <- NULL
  names(ddd_rates$mu_N) <- NULL
  # Test
  cat(paste("Testing ddmodel =", ddmodel, "\n"))
  expect_equal(ddd_rates, exptd_rates)
}

test_dd_lamuN <- function(ddmodel, pars_set, exptd_rates_set) {
  # global variables
  N <- pars_set["K"]
  x <- rep(NA, 10) # length of the Q_n vector
  n_seq <- c(0, 0:(length(x) + 2 * N)) # based on internal code, not sure why
  
  exptd_rates <- match_exptd_rates(
    ddmodel = ddmodel,
    exptd_rates_set = exptd_rates_set,
    n_seq = n_seq
  )
  ddd_rates <- list(
    "la_N" = unlist(lapply(n_seq, function(n) {
      dd_lamuN(ddmodel = ddmodel, pars = pars_set, N = n)[1]
    }), use.names = FALSE),
    "mu_N" = unlist(lapply(n_seq, function(n) {
      dd_lamuN(ddmodel = ddmodel, pars = pars_set, N = n)[2]
    }), use.names = FALSE)
  )
  # Test
  cat(paste("Testing ddmodel =", ddmodel, "\n"))
  expect_equal(ddd_rates, exptd_rates)
}

test_that("set1", {
  ddmodels <- c(5:8, 11:15)
  invisible(lapply(
    ddmodels, 
    test_dd_loglik_rhs_precomp, 
    pars_set = pars_set1, 
    exptd_rates_set = exptd_rates_set1
  ))
  invisible(lapply(
    ddmodels, 
    test_lambdamu, 
    pars_set = pars_set1, 
    exptd_rates_set = exptd_rates_set1
  ))
  invisible(lapply(
    ddmodels, 
    test_dd_lamuN, 
    pars_set = pars_set1, 
    exptd_rates_set = exptd_rates_set1
  ))
})

test_that("set2", {
  ddmodels <- c(1, 5:9, 11:15)
  invisible(lapply(
    ddmodels, 
    test_dd_loglik_rhs_precomp, 
    pars_set = pars_set2, 
    exptd_rates_set = exptd_rates_set2
  ))
  invisible(lapply(
    ddmodels, 
    test_lambdamu, 
    pars_set = pars_set2, 
    exptd_rates_set = exptd_rates_set2
  ))
  invisible(lapply(
    ddmodels, 
    test_dd_lamuN, 
    pars_set = pars_set2, 
    exptd_rates_set = exptd_rates_set2
  ))
})

test_that("set3", {
  ddmodels <- c(3, 5:8, 10:15)
  invisible(lapply(
    ddmodels, 
    test_dd_loglik_rhs_precomp, 
    pars_set = pars_set3, 
    exptd_rates_set = exptd_rates_set3
  ))
  invisible(lapply(
    ddmodels, 
    test_lambdamu, 
    pars_set = pars_set3, 
    exptd_rates_set = exptd_rates_set3
  ))
  invisible(lapply(
    ddmodels, 
    test_dd_lamuN, 
    pars_set = pars_set3, 
    exptd_rates_set = exptd_rates_set3
  ))
})

