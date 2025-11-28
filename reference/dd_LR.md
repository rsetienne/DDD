# Bootstrap likelihood ratio test of diversity-dependent diversification model

This function computes the maximum likelihood and the associated
estimates of the parameters of a diversity-dependent diversification
model for a given set of phylogenetic branching times. It then performs
a bootstrap likelihood ratio test of the diversity-dependent (DD) model
against the constant-rates (CR) birth-death model. Finally, it computes
the power of this test.

## Usage

``` r
dd_LR(
  brts,
  initparsoptDD,
  initparsoptCR,
  missnumspec,
  outputfilename = NULL,
  seed = 42,
  endmc = 1000,
  alpha = 0.05,
  plotit = TRUE,
  res = 10 * (1 + length(brts) + missnumspec),
  ddmodel = 1,
  cond = 1,
  btorph = 1,
  soc = 2,
  tol = c(0.001, 1e-04, 1e-06),
  maxiter = 2000,
  changeloglikifnoconv = FALSE,
  optimmethod = "simplex",
  methode = "analytical"
)
```

## Arguments

- brts:

  A set of branching times of a phylogeny, all positive

- initparsoptDD:

  The initial values of the parameters that must be optimized for the
  diversity-dependent (DD) model: lambda_0, mu and K

- initparsoptCR:

  The initial values of the parameters that must be optimized for the
  constant-rates (CR) model: lambda and mu

- missnumspec:

  The number of species that are in the clade but missing in the
  phylogeny

- outputfilename:

  The name (and location) of the file where the output will be saved.
  Default is no save.

- seed:

  The seed for the pseudo random number generator for simulating the
  bootstrap data

- endmc:

  The number of bootstraps

- alpha:

  The significance level of the test

- plotit:

  Boolean to plot results or not

- res:

  Sets the maximum number of species for which a probability must be
  computed, must be larger than 1 + length(brts)

- ddmodel:

  Sets the model of diversity-dependence:  
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
  `ddmodel == 5` : linear dependence in speciation and extinction rate  

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

- maxiter:

  Sets the maximum number of iterations in the optimization

- changeloglikifnoconv:

  if TRUE the loglik will be set to -Inf if ML does not converge

- optimmethod:

  Method used in optimization of the likelihood. Current default is
  'simplex'. Alternative is 'subplex' (default of previous versions)

- methode:

  The method used to solve the master equation, default is 'analytical'
  which uses matrix exponentiation; alternatively numerical ODE solvers
  can be used, such as 'odeint::runge_kutta_cash_karp54'. These were
  used in the package before version 3.1.

## Value

- treeCR:

  a list of trees generated under the constant-rates model using the ML
  parameters under the CR model

- treeDD:

  a list of trees generated under the diversity-dependent model using
  the ML parameters under the diversity-dependent model

- out:

  a dataframe with the parameter estimates and maximum likelihoods for
  diversity-dependent and constant-rates models `$model` - the model
  used to generate the data. 0 = unknown (for real data), 1 = CR, 2 =
  DD  
  `$mc` - the simulation number for each model  
  `$lambda_CR` - speciation rate estimated under CR  
  `$mu_CR` - extinction rate estimated under CR  
  `$LL_CR` - maximum likelihood estimated under CR  
  `$conv_CR` - convergence code for likelihood optimization; conv = 0
  means convergence  
  `$lambda_DD1` - initial speciation rate estimated under DD for first
  set of initial values  
  `$mu_DD1` - extinction rate estimated under DD for first set of
  initial values  
  `$K_DD1` - clade-wide carrying-capacity estimated under DD for first
  set of initial values  
  `$LL_DD1` - maximum likelihood estimated under DD for first set of
  initial values  
  `$conv_DD1` - convergence code for likelihood optimization for first
  set of initial values; conv = 0 means convergence  
  `$lambda_DD2` - initial speciation rate estimated under DD for second
  set of initial values  
  `$mu_DD2` - extinction rate estimated under DD for second set of
  initial values  
  `$K_DD2` - clade-wide carrying-capacity estimated under DD for second
  set of initial values  
  `$LL_DD2` - maximum likelihood estimated under DD for second set of
  initial values  
  `$conv_DD2` - convergence code for likelihood optimization for second
  set of initial values; conv = 0 means convergence  
  `$LR` - likelihood ratio between DD and CR

- pvalue:

  p-value of the test

- LRalpha:

  Likelihood ratio at the signifiance level alpha

- poweroftest:

  power of the test for significance level alpha

## Details

The output is a list with 3 elements:

## References

\- Etienne, R.S. et al. 2016. Meth. Ecol. Evol. 7: 1092-1099, doi:
10.1111/2041-210X.12565  
- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## See also

[`dd_loglik`](https://rsetienne.github.io/DDD/reference/dd_loglik.md),
[`dd_ML`](https://rsetienne.github.io/DDD/reference/dd_ML.md)

## Author

Rampal S. Etienne & Bart Haegeman
