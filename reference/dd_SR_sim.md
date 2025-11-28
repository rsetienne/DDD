# Function to simulate the diversity-dependent diversification process with a shift in one or more of the parameters

Simulating the diversity-dependent diversification process with a
parameter shift at a certain time

## Usage

``` r
dd_SR_sim(pars, age, ddmodel = 1)
```

## Arguments

- pars:

  Vector of parameters:  
    
  `pars[1]` corresponds to lambda1 (speciation rate before the rate
  shift)  
  `pars[2]` corresponds to mu1 (extinction rate before the rate shift)  
  `pars[3]` corresponds to K1 (clade-level carrying capacity before the
  rate shift)  
  `pars[4]` corresponds to lambda2 (speciation rate after the rate
  shift)  
  `pars[5]` corresponds to mu2 (extinction rate after the rate shift)  
  `pars[6]` corresponds to K2 (clade-level carrying capacity after the
  rate shift)  
  `pars[7]` corresponds to the time of shift

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
  - the fourth column is the time of extinction of the species If this
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
 dd_SR_sim(c(0.2,0.1,20,0.2,0.1,40,5),10) 
#> $tes
#> 
#> Phylogenetic tree with 3 tips and 2 internal nodes.
#> 
#> Tip labels:
#>   t5, t6, t7
#> 
#> Rooted; includes branch length(s).
#> 
#> $tas
#> 
#> Phylogenetic tree with 7 tips and 6 internal nodes.
#> 
#> Tip labels:
#>   t1, t5, t2, t3, t4, t6, ...
#> 
#> Rooted; includes branch length(s).
#> 
#> $L
#>            [,1] [,2] [,3]      [,4]
#> [1,] 10.0000000    0   -1  3.210579
#> [2,] 10.0000000   -1    2  5.349671
#> [3,]  8.8782487    2    3  3.928384
#> [4,]  6.7947575    3    4  1.747065
#> [5,]  6.3788456   -1   -5 -1.000000
#> [6,]  3.5658002    4    6 -1.000000
#> [7,]  0.1126442    6    7 -1.000000
#> 
#> $brts
#> [1] 10.0000000  0.1126442
#> 
```
