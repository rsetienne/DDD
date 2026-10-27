# Maximization of the loglikelihood under the diversity-independent, possibly time-dependent diversification model

This function computes the maximum likelihood estimates of the
parameters of a diversity-independent diversification model for a given
set of phylogenetic branching times. It also outputs the corresponding
loglikelihood that can be used in model comparisons.

## Usage

``` r
bd_ML(
  brts,
  initparsopt = c(0.1, 0.05 * (tdmodel <= 1) + 10 * (length(brts) + missnumspec) *
    (tdmodel > 1)),
  idparsopt = c(1, 2 + (tdmodel > 1)),
  idparsfix = (1:4)[-idparsopt],
  parsfix = rep(0, 4)[idparsfix],
  missnumspec = 0,
  tdmodel = 0,
  cond = 1,
  btorph = 1,
  soc = 2,
  tol = c(0.001, 1e-04, 1e-06),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  changeloglikifnoconv = FALSE,
  optimmethod = "simplex",
  num_cycles = 1,
  methode = "odeint::runge_kutta_cash_karp54",
  verbose = FALSE
)
```

## Arguments

- brts:

  A set of branching times of a phylogeny, all positive

- initparsopt:

  The initial values of the parameters that must be optimized

- idparsopt:

  The ids of the parameters that must be optimized, e.g. 1:3 for
  intrinsic speciation rate, extinction rate and carrying capacity. The
  ids are defined as follows:  
  id == 1 corresponds to lambda0 (speciation rate)  
  id == 2 corresponds to mu0 (extinction rate)  
  id == 3 corresponds to lamda1 (parameter controlling decline in
  speciation rate with time)  
  id == 4 corresponds to mu1 (parameter controlling decline in
  extinction rate with time)

- idparsfix:

  The ids of the parameters that should not be optimized, e.g. c(1,3) if
  lambda0 and lambda1 should not be optimized, but only mu0 and mu1. In
  that case idparsopt must be c(2,4). The default is to fix all
  parameters not specified in idparsopt.

- parsfix:

  The values of the parameters that should not be optimized

- missnumspec:

  The number of species that are in the clade but missing in the
  phylogeny

- tdmodel:

  Sets the model of time-dependence:  
  tdmodel == 0 : constant speciation and extinction rates  
  tdmodel == 1 : speciation and/or extinction rate is exponentially
  declining with time  
  tdmodel == 2 : stepwise decline in speciation rate as in
  diversity-dependence without extinction  
  tdmodel == 3 : decline in speciation rate following deterministic
  logistic equation for ddmodel = 1  
  tdmodel == 4 : decline in speciation rate such that the expected
  number of species matches with that of ddmodel = 1 with the same mu

- cond:

  Conditioning:  
  cond == 0 : conditioning on stem or crown age  
  cond == 1 : conditioning on stem or crown age and non-extinction of
  the phylogeny  
  cond == 2 : conditioning on stem or crown age and on the total number
  of extant taxa (including missing species)  
  cond == 3 : conditioning on the total number of extant taxa (including
  missing species)

- btorph:

  Sets whether the likelihood is for the branching times (0) or the
  phylogeny (1)

- soc:

  Sets whether stem or crown age should be used (1 or 2). Ignored when
  cond == 3: the likelihood is then computed by integrating over all
  possible stem ages, so the first branching time is always the crown
  age.

- tol:

  Sets the tolerances in the optimization. Consists of:  
  reltolx = relative tolerance of parameter values in optimization  
  reltolf = relative tolerance of function value in optimization  
  abstolx = absolute tolerance of parameter values in optimization

- maxiter:

  Sets the maximum number of iterations in the optimization

- changeloglikifnoconv:

  if TRUE the loglik will be set to -Inf if ML does not converge

- optimmethod:

  Method used in optimization of the likelihood. Current default is
  'simplex'. Alternative is 'subplex' (default of previous versions)

- num_cycles:

  the number of cycles of opimization. If set at Inf, it will do as many
  cycles as needed to meet the tolerance set for the target function.

- methode:

  The method used to solve the master equation under tdmodel = 4,
  default is 'odeint::runge_kutta_cash_karp54'.

- verbose:

  Sets the amount of output shown on screen:  
  verbose = 0 : no intermediate output (only the final result and error
  messages are shown)  
  verbose = 1 : which parameters are optimized and fixed, the
  loglikelihood of the initial parameter values and a trace of the
  iterations of the optimization  
  verbose = 2 : as 1, plus the parameters and loglikelihood for every
  call to the loglik function  
  verbose = 3 : as 2, plus all messages (e.g. about parameter values or
  numerical problems that lead to a loglikelihood of -Inf).  
  TRUE is interpreted as 1 and FALSE as 0.

## Value

A dataframe with the following elements:  

- lambda0:

  gives the maximum likelihood estimate of lambda0

- mu0:

  gives the maximum likelihood estimate of mu0

- lambda1:

  gives the maximum likelihood estimate of lambda1

- mu1:

  gives the maximum likelihood estimate of mu1

- loglik:

  gives the maximum loglikelihood

- df:

  gives the number of estimated parameters, i.e. degrees of feedom

- conv:

  gives a message on convergence of optimization; conv = 0 means
  convergence

## Details

The output is a dataframe containing estimated parameters and maximum
loglikelihood. The computed loglikelihood contains the factor q! m! /
(q + m)! where q is the number of species in the phylogeny and m is the
number of missing species, as explained in the supplementary material to
Etienne et al. 2012.

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## See also

[`bd_loglik`](https://rsetienne.github.io/DDD/reference/bd_loglik.md)

## Author

Rampal S. Etienne & Bart Haegeman

## Examples

``` r

cat("Estimating parameters for a set of branching times brts with the default settings:")
#> Estimating parameters for a set of branching times brts with the default settings:
brts = 1:20
bd_ML(brts = brts, cond = 1)
#> 
#> Maximum likelihood parameter estimates: lambda0: 0.082607, mu0: 0.000000, lambda1: 0.000000, mu1: 0.000000: 
#> Maximum loglikelihood: -66.379166
#>      lambda0          mu0 lambda1 mu1    loglik df conv
#> 1 0.08260722 5.813998e-16       0   0 -66.37917  2    0
```
