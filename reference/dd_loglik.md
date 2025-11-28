# Loglikelihood for diversity-dependent diversification models

This function computes loglikelihood of a diversity-dependent
diversification model for a given set of branching times and parameter
values.

## Usage

``` r
dd_loglik(pars1, pars2, brts, missnumspec, methode = "analytical")
```

## Arguments

- pars1:

  Vector of parameters:  
    
  `pars1[1]` corresponds to lambda (speciation rate)  
  `pars1[2]` corresponds to mu (extinction rate)  
  `pars1[3]` corresponds to K (clade-level carrying capacity)

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
  - `pars2[2] == 1.4` : positive diversity-dependence in speciation rate
  with parameter K' (= diversity where speciation rate reaches half its
  maximum); lambda = lambda0 \* S/(S + K') where S is species richness  
  - `pars2[2] == 1.5` : positive and negative diversity-dependence in
  speciation rate with parameter K' (= diversity where speciation = 0);
  lambda = lambda0 \* S/K' \* (1 - S/K') where S is species richness  
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
  - `pars2[2] == 5` linear dependence in speciation and extinction
  rate  
    
  `pars2[3]` sets the conditioning:  
  - `pars2[3] == 0` conditioning on stem or crown age  
  - `pars2[3] == 1` conditioning on stem or crown age and non-extinction
  of the phylogeny  
  - `pars2[3] == 2` conditioning on stem or crown age and on the total
  number of extant taxa (including missing species)  
  - `pars2[3] == 3` conditioning on the total number of extant taxa
  (including missing species)  
    
  `pars2[4]` sets whether the likelihood is for the branching times (0)
  or the phylogeny (1)  
    
  `pars2[5]` sets whether the parameters and likelihood should be shown
  on screen (1) or not (0)  
    
  `pars2[6]` sets whether the first data point is stem age (1) or crown
  age (2)  
    
  `pars2[7]` sets the absolute tolerance of the integration  
    
  `pars2[8]` sets the relative tolerance of the integration  
    
  `pars2[9]` sets threshold on the number of species for which
  logarithmic integration will be used. If the last three elements are
  not set ddefault values will be used of 1E-10, 1E-8 and Inf.

- brts:

  A set of branching times of a phylogeny, all positive

- missnumspec:

  The number of species that are in the clade but missing in the
  phylogeny

- methode:

  The method used to solve the master equation, default is 'analytical'
  which uses matrix exponentiation; alternatively numerical ODE solvers
  can be used. Before version 3.1 these were solvers from the deSolve
  package such as 'lsoda' and 'ode45'. Currently solvers from odeint are
  used, such as 'odeint::runge_kutta_cash_karp54',
  'odeint::runge_kutta_fehlberg78', 'odeint::runge_kutta_dopri5', or
  odeint::bulirsch_stoer'. The first two are recommended in most cases.

## Value

The loglikelihood

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## See also

[`dd_ML`](https://rsetienne.github.io/DDD/reference/dd_ML.md),
[`dd_SR_loglik`](https://rsetienne.github.io/DDD/reference/dd_SR_loglik.md),
[`dd_KI_loglik`](https://rsetienne.github.io/DDD/reference/dd_KI_loglik.md)

## Author

Rampal S. Etienne & Bart Haegeman

## Examples

``` r
dd_loglik(pars1 = c(0.5,0.1,100), pars2 = c(100,1,1,1,0,2), brts = 1:10, missnumspec = 0) 
#> [1] -34.13958
```
