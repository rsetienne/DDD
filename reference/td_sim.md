# Simulation of a diversity-dependent-like time-dependent process

Simulates a phylogenetic tree branching according to a time-dependent
process calibrated on the expected number of species under a
diversity-dependent process over time.

## Usage

``` r
td_sim(pars, age, ddmodel = 1, methode = "ode45")
```

## Arguments

- pars:

  Vector of parameters:  
    
  `pars[1]` corresponds to lambda0 (speciation rate)  
  `pars[2]` corresponds to mu0 (extinction rate)  
  `pars[3]` corresponds to lambda1 (decline parameter in speciation
  rate) or K in diversity-dependence-like models  
  `pars[4]` corresponds to mu1 (decline parameter in extinction rate)

- age:

  crown age of the tree to simulate, i.e. the simulation time.

- ddmodel:

  the diversity-dependent model used as reference for the time-dependent
  model.

- methode:

  The method used to solve the master equation. See
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html)
  documentation for possible inputs

## Value

A list with the following four elements: The first element is the tree
of extant species in phylo format  
The second element is the tree of all species, including extinct
species, in phylo format  
The third element is a matrix of all species where - the first column is
the time at which a species is born  
- the second column is the label of the parent of the species; positive
and negative values only indicate whether the species belongs to the
left or right crown lineage  
- the third column is the label of the daughter species itself; positive
and negative values only indicate whether the species belongs to the
left or right crown lineage  
- the fourth column is the time of extinction of the species. If this
equals -1, then the species is still extant.

## Author

César Martinez, Rampal S. Etienne
