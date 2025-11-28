# Loglikelihood for diversity-independent diversification model

This function computes loglikelihood of a diversity-independent
diversification model for a given set of branching times and parameter
values.

## Usage

``` r
bd_loglik(
  pars1,
  pars2,
  brts,
  missnumspec,
  methode = "odeint::runge_kutta_cash_karp54"
)
```

## Arguments

- pars1:

  Vector of parameters:  
    
  `pars1[1]` corresponds to lambda0 (speciation rate)  
  `pars1[2]` corresponds to mu0 (extinction rate)  
  `pars1[3]` corresponds to lambda1 (decline parameter in speciation
  rate) or K in diversity-dependence-like models  
  `pars1[4]` corresponds to mu1 (decline parameter in extinction rate)

- pars2:

  Vector of model settings:  
    
  `pars2[1]` sets the model of time-dependence:  
  - `pars2[1] == 0` no time dependence  
  - `pars2[1] == 1` speciation and/or extinction rate is exponentially
  declining with time  
  - `pars2[1] == 2` stepwise decline in speciation rate as in
  diversity-dependence without extinction  
  - `pars2[1] == 3` decline in speciation rate following deterministic
  logistic equation for ddmodel = 1  
  - `pars2[1] == 4` decline in speciation rate such that the expected
  number of species matches with that of ddmodel = 1 with the same mu  
    
  `pars2[2]` sets the conditioning:  
  - `pars[2] == 0` conditioning on stem or crown age  
  - `pars[2] == 1` conditioning on stem or crown age and non-extinction
  of the phylogeny  
  - `pars[2] == 2` conditioning on stem or crown age and on the total
  number of extant taxa (including missing species)  
  - `pars[2] == 3` conditioning on the total number of extant taxa
  (including missing species)  
    
  `pars2[3]` sets whether the likelihood is for the branching times (0)
  or the phylogeny (1)  
    
  `pars2[4]` sets whether the parameters and likelihood should be shown
  on screen (1) or not (0)  
    
  `pars2[5]` sets whether the first data point is stem age (1) or crown
  age (2)

- brts:

  A set of branching times of a phylogeny, all positive

- missnumspec:

  The number of species that are in the clade but missing in the
  phylogeny

- methode:

  The method used to solve the master equation, default is
  'odeint::runge_kutta_cash_karp54'.

## Value

The loglikelihood

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## See also

[`bd_ML`](https://rsetienne.github.io/DDD/reference/bd_ML.md)

## Author

Rampal S. Etienne, Bart Haegeman & Cesar Martinez

## Examples

``` r
bd_loglik(pars1 = c(0.5,0.1), pars2 = c(0,1,1,0,2), brts = 1:10, 
missnumspec = 0)
#> [1] -35.86764
```
