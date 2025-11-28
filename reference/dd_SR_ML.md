# Maximization of the loglikelihood under a diversity-dependent diversification model with a shift in the parameters

This function computes the maximum likelihood estimates of the
parameters of a diversity-dependent diversification model with shifting
parameters at time t = tshift for a given set of phylogenetic branching
times. It also outputs the corresponding loglikelihood that can be used
in model comparisons.

## Usage

``` r
dd_SR_ML(
  brts,
  initparsopt = c(0.5, 0.1, 2 * (1 + length(brts) + missnumspec), 2 * (1 + length(brts) +
    missnumspec), max(brts)/2),
  parsfix = NULL,
  idparsopt = c(1:3, 6:7),
  idparsfix = NULL,
  idparsnoshift = (1:7)[c(-idparsopt, (-1)^(length(idparsfix) != 0) * idparsfix)],
  res = 10 * (1 + length(brts) + missnumspec),
  ddmodel = 1,
  missnumspec = 0,
  cond = 1,
  btorph = 1,
  soc = 2,
  allbp = FALSE,
  tol = c(0.001, 1e-04, 1e-06),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  changeloglikifnoconv = FALSE,
  optimmethod = "simplex",
  num_cycles = 1,
  methode = "analytical",
  verbose = FALSE
)
```

## Arguments

- brts:

  A set of branching times of a phylogeny, all positive

- initparsopt:

  The initial values of the parameters that must be optimized

- parsfix:

  The values of the parameters that should not be optimized

- idparsopt:

  The ids of the parameters that must be optimized, e.g. 1:7 for all
  parameters. The ids are defined as follows:  
  id == 1 corresponds to lambda (speciation rate) before the shift  
  id == 2 corresponds to mu (extinction rate) before the shift  
  id == 3 corresponds to K (clade-level carrying capacity) before the
  shift  
  id == 4 corresponds to lambda (speciation rate) after the shift  
  id == 5 corresponds to mu (extinction rate) after the shift  
  id == 6 corresponds to K (clade-level carrying capacity) after the
  shift  
  id == 7 corresponds to tshift (the time of shift)

- idparsfix:

  The ids of the parameters that should not be optimized, e.g.
  c(1,3,4,6) if lambda and K should not be optimized, but only mu. In
  that case idparsopt must be c(2,5,7). The default is to fix all
  parameters not specified in idparsopt.

- idparsnoshift:

  The ids of the parameters that should not shift; This can only apply
  to ids 4, 5 and 6, e.g. idparsnoshift = c(4,5) means that lambda and
  mu have the same values before and after tshift

- res:

  sets the maximum number of species for which a probability must be
  computed, must be larger than 1 + length(brts)

- ddmodel:

  sets the model of diversity-dependence:  
  ddmodel == 1 : linear dependence in speciation rate  
  ddmodel == 2 : exponential dependence in speciation rate  
  ddmodel == 2.1 : variant of exponential dependence in speciation rate
  with offset at infinity  
  ddmodel == 2.2 : 1/n dependence in speciation rate  
  ddmodel == 3 : linear dependence in extinction rate  
  ddmodel == 4 : exponential dependence in extinction rate  
  ddmodel == 4.1 : variant of exponential dependence in extinction rate
  with offset at infinity  
  ddmodel == 4.2 : 1/n dependence in extinction rate with offset at
  infinity

- missnumspec:

  The number of species that are in the clade but missing in the
  phylogeny

- cond:

  Conditioning:  
  cond == 0 : no conditioning  
  cond == 1 : conditioning on non-extinction of the phylogeny  
  cond == 2 : conditioning on non-extinction of the phylogeny and on the
  total number of extant taxa (including missing species)  
  cond == 3 : conditioning on the total number of extant taxa (including
  missing species)  
    
  Note: cond == 3 assumes a uniform prior on stem age, as is the
  standard in constant-rate birth-death models, see e.g. D. Aldous & L.
  Popovic 2004. Adv. Appl. Prob. 37: 1094-1115 and T. Stadler 2009. J.
  Theor. Biol. 261: 58-66.

- btorph:

  Sets whether the likelihood is for the branching times (0) or the
  phylogeny (1)

- soc:

  Sets whether stem or crown age should be used (1 or 2)

- allbp:

  Sets whether a search should be done with various initial conditions,
  with tshift at each of the branching points (TRUE/FALSE)

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

  The method used to solve the master equation, default is 'analytical'
  which uses matrix exponentiation; alternatively numerical ODE solvers
  can be used, such as 'odeint::runge_kutta_cash_karp54'. These were
  used in the package before version 3.1.

- verbose:

  Show the parameters and loglikelihood for every call to the loglik
  function

## Value

- lambda_1:

  gives the maximum likelihood estimate of lambda before the shift

- mu_1:

  gives the maximum likelihood estimate of mu before the shift

- K_1:

  gives the maximum likelihood estimate of K before the shift

- lambda_2:

  gives the maximum likelihood estimate of lambda after the shift

- mu_2:

  gives the maximum likelihood estimate of mu after the shift

- K_2:

  gives the maximum likelihood estimate of K after the shift

- t_shift:

  gives the time of the shift

- loglik:

  gives the maximum loglikelihood

- df:

  gives the number of estimated parameters, i.e. degrees of feedom

- conv:

  gives a message on convergence of optimization; conv = 0 means
  convergence

## Details

The output is a dataframe containing estimated parameters and maximum
loglikelihood. The computed loglikelihood contains the factor q! m!/(q +
m)! where q is the number of species in the phylogeny and m is the
number of missing species, as explained in the supplementary material to
Etienne et al. 2012.

## Note

The optimization may get trapped in local optima. Try different starting
values to search for the global optimum.

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## See also

[`dd_SR_loglik`](https://rsetienne.github.io/DDD/reference/dd_SR_loglik.md),
[`dd_ML`](https://rsetienne.github.io/DDD/reference/dd_ML.md),
[`dd_KI_ML`](https://rsetienne.github.io/DDD/reference/dd_KI_ML.md),

## Author

Rampal S. Etienne & Bart Haegeman

## Examples

``` r
cat("This will estimate parameters for a sets of branching times brts without conditioning.\n")
#> This will estimate parameters for a sets of branching times brts without conditioning.
cat("The tolerance of the optimization is set ridiculously high to make runtime fast.\n")
#> The tolerance of the optimization is set ridiculously high to make runtime fast.
cat("In real applications, use the default or more stringent settings for tol.\n")
#> In real applications, use the default or more stringent settings for tol.
brts = 1:10
dd_SR_ML(brts = brts, initparsopt = c(0.4581, 1E-6, 17.69, 11.09, 8.9999), idparsopt = c(1:3,6,7),
         idparsfix = NULL, parsfix = NULL, idparsnoshift = c(4,5), cond = 0,
         tol = c(1E-1,1E-1,1E-1),optimmethod = 'simplex'
)
#> You are optimizing la mu K K2 tshift 
#> You are fixing nothing 
#> You are not shifting la2 mu2 
#> Optimizing the likelihood - this may take a while. 
#> The loglikelihood for the initial parameter values is -24.52893 
#> 1 0.4581 1e-06 17.69 11.09 8.9999 -24.5289252582909 initial 
#> Optimization has terminated successfully. 
#> 
#> Maximum likelihood parameter estimates: 0.458100 0.000001 18.574500 0.458100 0.000001 11.090000 8.999900
#> Maximum loglikelihood: -24.527807
#> 
```
