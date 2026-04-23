# Function to simulate a key innovation in macro-evolution with the innovative clade decoupling from the diversity-dependent diversification dynamics of the main clade

Simulating a diversity-dependent diversification process where at a
given time a new clade emerges with different inherent speciation rate
and extinction rate and clade-level carrying capacity and with decoupled
dynamics

## Usage

``` r
dd_KI_sim(pars, age, ddmodel = 1)
```

## Arguments

- pars:

  Vector of parameters:  
    
  `pars[1]` corresponds to lambda_M (speciation rate of the main
  clade)  
  `pars[2]` corresponds to mu_M (extinction rate of the main clade)  
  `pars[3]` corresponds to K_M (clade-level carrying capacity of the
  main clade) `pars[4]` corresponds to lambda_S (speciation rate of the
  subclade)  
  `pars[5]` corresponds to mu_S (extinction rate of the subclade)  
  `pars[5]` corresponds to K_S (clade-level carrying capacity of the
  subclade)  
  `pars[7]` tinn, the time the shift in rates occurs in the lineage
  leading to the subclade

- age:

  Sets the crown age for the simulation

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

## Value

- out :

  A list with the following elements: The first element is the tree of
  extant species in phylo format  
  The second element is the tree of all species, including extinct
  species, in phylo format  
  The third element is a matrix of all species where  
  - the first column is the time at which a species is born  
  - the second column is the label of the parent of the species;
  positive and negative values only indicate whether the species belongs
  to the left or right crown lineage  
  - the third column is the label of the daughter species itself;
  positive and negative values only indicate whether the species belongs
  to the left or right crown lineage  
  - the fourth column is the time of extinction of the species  
  If the fourth element equals -1, then the species is still extant.  
  - the fifth column indicates whether the species belong to the main
  clade (0) or the subclade (1)  
  The fourth element is the subclade tree of extant species (without
  stem)  
  The fifth element is the subclade tree of all species (without stem)  
  The sixth element is the same as the first, except that it has
  attributed 0 for the main clade and 1 for the subclade  
  The seventh element is the same as the Second, except that it has
  attributed 0 for the main clade and 1 for the subclade  
  The sixth and seventh element will be NULL if the subclade does not
  exist (because it went extinct).

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## Author

Rampal S. Etienne

## Examples

``` r
 dd_KI_sim(c(0.2,0.1,20,0.1,0.05,30,4),10) 


#> $tes
#> 
#> Phylogenetic tree with 6 tips and 5 internal nodes.
#> 
#> Tip labels:
#>   t1, t7, t9, t8, t3, t4
#> 
#> Rooted; includes branch length(s).
#> 
#> $tas
#> 
#> Phylogenetic tree with 9 tips and 8 internal nodes.
#> 
#> Tip labels:
#>   t1, t7, t9, t8, t5, t2, ...
#> 
#> Rooted; includes branch length(s).
#> 
#> $L
#>            [,1] [,2] [,3]      [,4] [,5]
#>  [1,] 10.000000    0   -1 -1.000000    0
#>  [2,] 10.000000   -1    2  8.003286    0
#>  [3,]  9.381329    2    3 -1.000000    0
#>  [4,]  7.972618    3    4 -1.000000    1
#>  [5,]  7.650547   -1   -5  5.786228    0
#>  [6,]  6.857887    3    6  2.240755    0
#>  [7,]  2.054015   -1   -7 -1.000000    0
#>  [8,]  1.426913   -7   -8 -1.000000    0
#>  [9,]  1.011446   -7   -9 -1.000000    0
#> 
#> $tesS
#> NULL
#> 
#> $tasS
#> NULL
#> 
#> $tes2
#> 
#> Phylogenetic tree with 6 tips and 5 internal nodes.
#> 
#> Tip labels:
#>  t1, t7, t9, t8, t3, t4
#> The tree includes a mapped, 2-state discrete character
#> with states:
#>  0, 1
#> 
#> Rooted; includes branch lengths.
#> 
#> $tas2
#> 
#> Phylogenetic tree with 9 tips and 8 internal nodes.
#> 
#> Tip labels:
#>  t1, t7, t9, t8, t5, t2, ...
#> 
#> The tree includes a mapped, 2-state discrete character
#> with states:
#>  0, 1
#> 
#> Rooted; includes branch lengths.
#> 
```
