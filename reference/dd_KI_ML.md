# Maximization of the loglikelihood under a diversity-dependent diversification model with decoupling of a subclade's diversication dynamics from the main clade's dynamics

This function computes the maximum likelihood estimates of the
parameters of a diversity-dependent diversification model with
decoupling of the diversification dynamics of a subclade from the
dynamics of the main clade for a given set of phylogenetic branching
times of main clade and subclade and the time of splitting of the
lineage that will form the subclade. It also outputs the corresponding
loglikelihood that can be used in model comparisons.

## Usage

``` r
dd_KI_ML(
  brtsM,
  brtsS,
  tsplit,
  initparsopt = c(0.5, 0.1, 2 * (1 + length(brtsM) + missnumspec[1]), 2 * (1 +
    length(brtsS) + missnumspec[length(missnumspec)]), (tsplit + max(brtsS))/2),
  parsfix = NULL,
  idparsopt = c(1:3, 6:7),
  idparsfix = NULL,
  idparsnoshift = (1:7)[c(-idparsopt, (-1)^(length(idparsfix) != 0) * idparsfix)],
  res = 10 * (1 + length(c(brtsM, brtsS)) + sum(missnumspec)),
  ddmodel = 1,
  missnumspec = 0,
  cond = 1,
  soc = 2,
  tol = c(0.001, 1e-04, 1e-06),
  maxiter = 1000 * round((1.25)^length(idparsopt)),
  changeloglikifnoconv = FALSE,
  optimmethod = "simplex",
  num_cycles = 1,
  methode = "analytical",
  correction = TRUE,
  verbose = FALSE
)
```

## Arguments

- brtsM:

  A set of branching times of the main clade in a phylogeny, all
  positive

- brtsS:

  A set of branching times of the subclade in a phylogeny, all positive

- tsplit:

  The branching time at which the lineage forming the subclade branches
  off, positive

- initparsopt:

  The initial values of the parameters that must be optimized

- parsfix:

  The values of the parameters that should not be optimized

- idparsopt:

  The ids of the parameters that must be optimized, e.g. 1:7 for all
  parameters. The ids are defined as follows:  
  id == 1 corresponds to lambda_M (speciation rate) of the main clade  
  id == 2 corresponds to mu_M (extinction rate) of the main clade  
  id == 3 corresponds to K_M (clade-level carrying capacity) of the main
  clade  
  id == 4 corresponds to lambda_S (speciation rate) of the subclade  
  id == 5 corresponds to mu_S (extinction rate) of the subclade  
  id == 6 corresponds to K_S (clade-level carrying capacity) of the
  subclade  
  id == 7 corresponds to t_d (the time of decoupling)

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
  computed, must be larger than 1 + max(length(brtsM),length(brtsS))

- ddmodel:

  sets the model of diversity-dependence:  
  `ddmodel == 1` : linear dependence in speciation rate with parameter K
  (= diversity where speciation = extinction)  
  `ddmodel == 1.3` : linear dependence in speciation rate with parameter
  K' (= diversity where speciation = 0)  
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

- missnumspec:

  The number of species that are in the clade but missing in the
  phylogeny. One can specify the sum of the missing species in main
  clade and subclade or a vector c(missnumspec_M,missnumspec_S) with
  missing species in main clade and subclade respectively.

- cond:

  Conditioning:  
  cond == 0 : no conditioning  
  cond == 1 : conditioning on non-extinction of the phylogeny  

- soc:

  Sets whether stem or crown age should be used (1 or 2); stem age only
  works when cond = 0

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

- correction:

  Sets whether the correction should be applied (TRUE) or not (FALSE)

- verbose:

  Show the parameters and loglikelihood for every call to the loglik
  function

## Value

- lambda_M:

  gives the maximum likelihood estimate of lambda of the main clade

- mu_M:

  gives the maximum likelihood estimate of mu of the main clade

- K_M:

  gives the maximum likelihood estimate of K of the main clade

- lambda_2:

  gives the maximum likelihood estimate of lambda of the subclade

- mu_S:

  gives the maximum likelihood estimate of mu of the subclade

- K_S:

  gives the maximum likelihood estimate of K of the subclade

- t_d:

  gives the time of the decoupling event

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

[`dd_KI_loglik`](https://rsetienne.github.io/DDD/reference/dd_KI_loglik.md),
[`dd_ML`](https://rsetienne.github.io/DDD/reference/dd_ML.md),
[`dd_SR_ML`](https://rsetienne.github.io/DDD/reference/dd_SR_ML.md),

## Author

Rampal S. Etienne & Bart Haegeman

## Examples

``` r
cat("This will estimate parameters for two sets of branching times brtsM, brtsS\n")
#> This will estimate parameters for two sets of branching times brtsM, brtsS
cat("without conditioning.\n")
#> without conditioning.
cat("The tolerance of the optimization is set high so runtime is fast in this example.\n")
#> The tolerance of the optimization is set high so runtime is fast in this example.
cat("In real applications, use the default or more stringent settins for tol.\n")
#> In real applications, use the default or more stringent settins for tol.
brtsM = 4:10
brtsS = seq(0.1,3.5,0.7)
tsplit = 5
dd_KI_ML(brtsM = brtsM, brtsS = brtsS, tsplit = tsplit, idparsopt = c(1:3,6,7),
  initparsopt = c(0.885, 2e-14, 6.999, 6.848, 4.001), idparsfix = NULL,
  parsfix = NULL,idparsnoshift = c(4,5), cond = 0, tol = c(3E-1,3E-1,3E-1),
  optimmethod = 'simplex')
#> You are optimizing la_M mu_M K_M K_S t_d 
#> You are fixing nothing 
#> You are not shifting la_S mu_S 
#> Optimizing the likelihood - this may take a while. 
#> The loglikelihood for the initial parameter values is -24.45333 
#> 1 0.885 2e-14 6.999 6.848 4.001 -24.4533303187317 initial 
#> Optimization has terminated successfully. 
#> 
#> Maximum likelihood parameter estimates: 0.885000 0.000000 6.999000 0.885000 0.000000 6.848000 4.001000
#> Maximum loglikelihood: -24.453330
```
