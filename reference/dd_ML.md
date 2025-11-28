# Maximization of the loglikelihood under a diversity-dependent diversification model

This function computes the maximum likelihood estimates of the
parameters of a diversity-dependent diversification model for a given
set of phylogenetic branching times. It also outputs the corresponding
loglikelihood that can be used in model comparisons.

## Usage

``` r
dd_ML(
  brts,
  initparsopt = initparsoptdefault(ddmodel, brts, missnumspec),
  idparsopt = 1:length(initparsopt),
  idparsfix = (1:(3 + (ddmodel == 5)))[-idparsopt],
  parsfix = parsfixdefault(ddmodel, brts, missnumspec, idparsopt),
  res = 10 * (1 + length(brts) + missnumspec),
  ddmodel = 1,
  missnumspec = 0,
  cond = 1,
  btorph = 1,
  soc = 2,
  tol = c(0.001, 1e-04, 1e-06),
  tolint = c(1e-10, 1e-08),
  probs_threshold = 0,
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

- idparsopt:

  The ids of the parameters that must be optimized, e.g. 1:3 for
  intrinsic speciation rate, extinction rate and carrying capacity. The
  ids are defined as follows:  
  id == 1 corresponds to lambda (speciation rate)  
  id == 2 corresponds to mu (extinction rate)  
  id == 3 corresponds to K (clade-level carrying capacity)  
  id == 4 corresponds to r (r = b/a where mu = mu_0 + b \* N and lambda
  = lambda_0 - a \* N) (This is only available when ddmodel = 5)

- idparsfix:

  The ids of the parameters that should not be optimized, e.g. c(1,3) if
  lambda and K should not be optimized, but only mu. In that case
  idparsopt must be 2. The default is to fix all parameters not
  specified in idparsopt.

- parsfix:

  The values of the parameters that should not be optimized

- res:

  Sets the maximum number of species for which a probability must be
  computed, must be larger than 1 + length(brts)

- ddmodel:

  Sets the model of diversity-dependence:  
  `ddmodel == 1` : linear dependence in speciation rate with parameter K
  (= diversity where speciation = extinction)  
  `ddmodel == 1.3` : linear dependence in speciation rate with parameter
  K' (= diversity where speciation = 0)  
  `ddmodel == 1.4` : positive diversity-dependence in speciation rate
  with parameter K' (= diversity where speciation rate reaches half its
  maximum); lambda = lambda0 \* S/(S + K') where S is species richness  
  `ddmodel == 1.5` : positive and negative dependence in speciation rate
  with parameter K' (= diversity where speciation = 0); lambda = lambda0
  \* S/K' \* (1 - S/K') where S is species richness  
  `ddmodel == 2` : exponential dependence in speciation rate with
  parameter K (= diversity where speciation = extinction)  
  `ddmodel == 2.1` : variant of exponential dependence in speciation
  rate with offset at infinity  
  `ddmodel == 2.2` : 1/n dependence in speciation rate  
  `ddmodel == 2.3` : exponential dependence in speciation rate with
  parameter x (= exponent)  
  `ddmodel == 3` : linear dependence in extinction rate  
  `ddmodel == 4` : exponential dependence in extinction rate  
  `ddmodel == 4.1` : variant of exponential dependence in extinction
  rate with offset at infinity  
  `ddmodel == 4.2` : 1/n dependence in extinction rate with offset at
  infinity  
  `ddmodel == 5` : linear dependence in speciation and extinction rate  

- missnumspec:

  The number of species that are in the clade but missing in the
  phylogeny

- cond:

  Conditioning:  
  cond == 0 : conditioning on stem or crown age  
  cond == 1 : conditioning on stem or crown age and non-extinction of
  the phylogeny  
  cond == 2 : conditioning on stem or crown age and on the total number
  of extant taxa (including missing species)  
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

- tol:

  Sets the tolerances in the optimization. Consists of:  
  reltolx = relative tolerance of parameter values in optimization  
  reltolf = relative tolerance of function value in optimization  
  abstolx = absolute tolerance of parameter values in optimization

- tolint:

  Sets the tolerance of the numerical integration. Consists of:  
  absoltint = absolute tolerance and  
  reltolint = relative tolerance.

- probs_threshold:

  Sets the threshold of the probability below which logarithmic
  integration must be used. Default is 0.

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

- lambda:

  gives the maximum likelihood estimate of lambda

- mu:

  gives the maximum likelihood estimate of mu

- K:

  gives the maximum likelihood estimate of K

- r:

  (only if ddmodel == 5) gives the ratio of linear dependencies in
  speciation and extinction rates

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

[`dd_loglik`](https://rsetienne.github.io/DDD/reference/dd_loglik.md),
[`dd_SR_ML`](https://rsetienne.github.io/DDD/reference/dd_SR_ML.md),
[`dd_KI_ML`](https://rsetienne.github.io/DDD/reference/dd_KI_ML.md),

## Author

Rampal S. Etienne & Bart Haegeman

## Examples

``` r
cat("Estimating the intrinsic speciation rate lambda and the carrying capacity K")
#> Estimating the intrinsic speciation rate lambda and the carrying capacity K
cat("for a fixed extinction rate of 0.1, conditioning on clade survival and two missing species:")
#> for a fixed extinction rate of 0.1, conditioning on clade survival and two missing species:
brts = 1:5
dd_ML(brts = brts,initparsopt = c(1.3078,7.4188), idparsopt = c(1,3), parsfix = 0.1,
      cond = 1, missnumspec = 2, tol = c(1E-3,1E-3,1E-4), optimmethod = 'simplex')
#> You are optimizing lambda K 
#> You are fixing mu 
#> Optimizing the likelihood - this may take a while. 
#> The loglikelihood for the initial parameter values is -9.349088 
#> 1 1.3078 7.41879999999999 -9.3490879225759 initial 
#> 2 1.3078 7.41879999999999 -9.3490879225759 contract inside 
#> 3 1.3078 7.41879999999999 -9.3490879225759 contract inside 
#> 4 1.3078 7.41879999999999 -9.3490879225759 contract inside 
#> 5 1.3078 7.41879999999999 -9.3490879225759 reflect 
#> 6 1.3078 7.41879999999999 -9.3490879225759 contract inside 
#> 7 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 8 1.3078 7.41879999999999 -9.3490879225759 contract inside 
#> 9 1.3078 7.41879999999999 -9.3490879225759 reflect 
#> 10 1.3078 7.41879999999999 -9.3490879225759 contract inside 
#> 11 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 12 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 13 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 14 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 15 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 16 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 17 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 18 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 19 1.3078 7.41879999999999 -9.3490879225759 contract outside 
#> 20 1.3078 7.41879999999999 -9.3490879225759 contract inside 
#> Optimization has terminated successfully. 
#> 
#> Maximum likelihood parameter estimates: lambda: 1.307800, mu: 0.100000, K: 7.418800
#> Maximum loglikelihood: -9.349088
```
