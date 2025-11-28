# Function to simulate the diversity-dependent diversification process

Simulating the diversity-dependent diversification process

## Usage

``` r
dd_sim(pars, age, ddmodel = 1)
```

## Arguments

- pars:

  Vector of parameters:  
    
  `pars[1]` corresponds to lambda (speciation rate)  
  `pars[2]` corresponds to mu (extinction rate)  
  `pars[3]` corresponds to K (clade-level carrying capacity)

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
  `ddmodel == 5` : linear dependence in speciation and extinction rate

## Value

- out :

  A list with the following four elements: The first element is the tree
  of extant species in phylo format  
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
  - the fourth column is the time of extinction of the species. If this
  equals -1, then the species is still extant.  
  The fourth element is the set of branching times of the tree of extant
  species.  

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## Author

Rampal S. Etienne

## Examples

``` r
 dd_sim(c(0.2,0.1,20),10) 
#> $tes
#> 
#> Phylogenetic tree with 2 tips and 1 internal node.
#> 
#> Tip labels:
#>   t1, t2
#> 
#> Rooted; includes branch length(s).
#> 
#> $tas
#> 
#> Phylogenetic tree with 3 tips and 2 internal nodes.
#> 
#> Tip labels:
#>   t1, t3, t2
#> 
#> Rooted; includes branch length(s).
#> 
#> $L
#>           [,1] [,2] [,3]      [,4]
#> [1,] 10.000000    0   -1 -1.000000
#> [2,] 10.000000   -1    2 -1.000000
#> [3,]  1.584906   -1   -3  1.241685
#> 
#> $brts
#> [1] 10
#> 
```
