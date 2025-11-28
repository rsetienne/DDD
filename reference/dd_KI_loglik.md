# Loglikelihood for diversity-dependent diversification models with decoupling of a subclade from a main clade at time t = t_d

This function computes loglikelihood of a diversity-dependent
diversification model for a given set of branching times and parameter
values where the diversity-dependent dynamics of a subclade decouple
from the dynamics of the main clade at time t_d, potentially accompanied
by a shift in parameters.

## Usage

``` r
dd_KI_loglik(
  pars1,
  pars2,
  brtsM,
  brtsS,
  missnumspec,
  methode = "odeint::runge_kutta_cash_karp54"
)
```

## Arguments

- pars1:

  Vector of parameters:  
    
  `pars1[1]` corresponds to lambda_M (speciation rate) of the main
  clade  
  `pars1[2]` corresponds to mu_M (extinction rate) of the main clade  
  `pars1[3]` corresponds to K_M (clade-level carrying capacity) of the
  main clade  
  `pars1[4]` corresponds to lambda_S (speciation rate) of the subclade  
  `pars1[5]` corresponds to mu_S (extinction rate) of the subclade  
  `pars1[6]` corresponds to K_S (clade-level carrying capacity) of the
  subclade  
  `pars1[7]` corresponds to t_d (the time of decoupling)

- pars2:

  Vector of model settings:  
    
  `pars2[1]` sets the maximum number of species for which a probability
  must be computed. This must be larger than 1 + missnumspec +
  length(brts).  
    
  `pars2[2]` sets the model of diversity-dependence:  
  - `pars2[2] == 1` linear dependence in speciation rate with parameter
  K (= diversity where speciation = extinction)  
  - `pars2[2] == 1.3` linear dependence in speciation rate with
  parameter K' (= diversity where speciation = 0)  
  - `pars2[2] == 2` exponential dependence in speciation rate with
  parameter K (= diversity where speciation = extinction)  
  - `pars2[2] == 2.1` variant of exponential dependence in speciation
  rate with offset at infinity  
  - `pars2[2] == 2.2` 1/n dependence in speciation rate  
  - `pars2[2] == 2.3` exponential dependence in speciation rate with
  parameter x (= exponent)  
  - `pars2[2] == 3` linear dependence in extinction rate  
  - `pars2[2] == 4` exponential dependence in extinction rate  
  - `pars2[2] == 4.1` variant of exponential dependence in extinction
  rate with offset at infinity  
  - `pars2[2] == 4.2` 1/n dependence in extinction rate  
    
  `pars2[3]` sets the conditioning:  
  - `pars2[3] == 0` no conditioning (or just crown age)  
  - `pars2[3] == 1` conditioning on non-extinction of the phylogeny  
  - `pars2[3] == 2` conditioning on number of species and crown age; not
  yet implemented  
  - `pars2[3] == 3` conditioning on number of species only; not yet
  implemented  
  - `pars2[3] == 4` conditioning on survival of the subclade  
  - `pars2[3] == 5` conditioning on survival of all subclades and of
  both crown lineages in the main clade. This assumes that subclades
  that have already shifted do not undergo another shift, i.e. shifts
  only occur in the main clade.  
    
  `pars2[4]` Obsolete.  
    
  `pars2[5]` sets whether the parameters and likelihood should be shown
  on screen (1) or not (0)  
    
  `pars2[6]` sets whether the first data point is stem age (1) or crown
  age (2)  
    
  `pars2[7]` sets whether the old (incorrect) likelihood should be used
  (0), or whether the new corrected likelihood should be used (1).

- brtsM:

  A set of branching times of the main clade in the phylogeny, all
  positive

- brtsS:

  A set of branching times of the subclade in the phylogeny, all
  positive

- missnumspec:

  The number of species that are in the clade but missing in the
  phylogeny. One can specify the sum of the missing species in main
  clade and subclade or a vector c(missnumspec_M,missnumspec_S) with
  missing species in main clade and subclade respectively.

- methode:

  The method used to solve the master equation, default is 'analytical'
  which uses matrix exponentiation; alternatively numerical ODE solvers
  can be used, such as 'odeint::runge_kutta_cash_karp54'. These were
  used in the package before version 3.1.

## Value

The loglikelihood

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## See also

[`dd_KI_ML`](https://rsetienne.github.io/DDD/reference/dd_KI_ML.md),
[`dd_loglik`](https://rsetienne.github.io/DDD/reference/dd_loglik.md)
[`dd_SR_loglik`](https://rsetienne.github.io/DDD/reference/dd_SR_loglik.md)

## Author

Rampal S. Etienne & Bart Haegeman

## Examples

``` r
pars1 = c(0.25,0.12,25.51,1.0,0.16,8.61,9.8)
pars2 = c(200,1,0,18.8,1,2)
missnumspec = 0
brtsM = c(25.2,24.6,24.0,22.5,21.7,20.4,19.9,19.7,18.8,17.1,15.8,11.8,9.7,8.9,5.7,5.2)
brtsS = c(9.6,8.6,7.4,4.9,2.5)
dd_KI_loglik(pars1,pars2,brtsM,brtsS,missnumspec)
#> Parameters: 0.250000 0.120000 25.510000 1.000000 0.160000 8.610000 9.800000, Loglikelihood: -77.815687
#> [1] -77.81569
```
