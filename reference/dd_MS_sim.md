# Function to simulate the macro-evolutionary succession process assuming diversity-dependent diversification

Simulating a diversity-dependent diversification process where at a
given time a new clade emerges with different inherent speciation rate
and extinction rate

## Usage

``` r
dd_MS_sim(pars, age, ddmodel = 1.3)
```

## Arguments

- pars:

  Vector of parameters:  
    
  `pars[1]` corresponds to lambda_M (speciation rate of the main
  clade)  
  `pars[2]` corresponds to mu_M (extinction rate of the main clade)  
  `pars[3]` corresponds to K' (maximum number of species or a proxy for
  it in case of exponential decline in speciation rate) `pars[4]`
  corresponds to lambda_S (speciation rate of the novel subclade)  
  `pars[5]` corresponds to mu_S (extinction rate)  
  `pars[6]` tinn, the time the shift in rates occurs in the lineage
  leading to the subclade

- age:

  Sets the crown age for the simulation

- ddmodel:

  Sets the model of diversity-dependence:  
  `ddmodel == 1.3` : linear dependence in speciation rate with parameter
  K' (= diversity where speciation = 0); ddmodel = 1 will be interpreted
  as this model  
  `ddmodel == 2.1` : variant of exponential dependence in speciation
  rate with offset at infinity; ddmodel = 2 will be interpreted as this
  model  
  `ddmodel == 2.2` : 1/n dependence in speciation rate  
  `ddmodel == 2.3` : exponential dependence in speciation rate with
  parameter x (= exponent)

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
 dd_MS_sim(c(0.2,0.1,20,0.1,0.05,4),10) 


#> $tes
#> 
#> Phylogenetic tree with 7 tips and 6 internal nodes.
#> 
#> Tip labels:
#>   t3, t7, t2, t10, t6, t8, ...
#> 
#> Rooted; includes branch length(s).
#> 
#> $tas
#> 
#> Phylogenetic tree with 11 tips and 10 internal nodes.
#> 
#> Tip labels:
#>   t1, t3, t7, t2, t10, t6, ...
#> 
#> Rooted; includes branch length(s).
#> 
#> $L
#>            [,1] [,2] [,3]       [,4] [,5]
#>  [1,] 10.000000    0   -1  8.2179313    0
#>  [2,] 10.000000   -1    2 -1.0000000    0
#>  [3,]  9.937863   -1   -3 -1.0000000    1
#>  [4,]  7.296311    2    4  3.5900766    0
#>  [5,]  4.068823    4    5  3.6695738    0
#>  [6,]  3.340208    2    6 -1.0000000    0
#>  [7,]  3.156793   -3   -7 -1.0000000    1
#>  [8,]  3.124882    6    8 -1.0000000    0
#>  [9,]  1.808820    6    9  0.9247604    0
#> [10,]  1.683834    2   10 -1.0000000    0
#> [11,]  1.209194    8   11 -1.0000000    0
#> 
#> $tesS
#> 
#> Phylogenetic tree with 2 tips and 1 internal node.
#> 
#> Tip labels:
#>   t3, t7
#> 
#> Rooted; includes branch length(s).
#> 
#> $tasS
#> 
#> Phylogenetic tree with 2 tips and 1 internal node.
#> 
#> Tip labels:
#>   t3, t7
#> 
#> Rooted; includes branch length(s).
#> 
#> $tes2
#> 
#> Phylogenetic tree with 7 tips and 6 internal nodes.
#> 
#> Tip labels:
#>  t3, t7, t2, t10, t6, t8, ...
#> 
#> The tree includes a mapped, 2-state discrete character
#> with states:
#>  0, 1
#> 
#> Rooted; includes branch lengths.
#> 
#> $tas2
#> 
#> Phylogenetic tree with 11 tips and 10 internal nodes.
#> 
#> Tip labels:
#>  t1, t3, t7, t2, t10, t6, ...
#> 
#> The tree includes a mapped, 2-state discrete character
#> with states:
#>  0, 1
#> 
#> Rooted; includes branch lengths.
#> 
```
