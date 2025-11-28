# Function to convert a table with speciation and extinction events to a phylogeny

Converting a table with speciation and extinction events to a phylogeny

## Usage

``` r
L2phylo(L, dropextinct = T)
```

## Arguments

- L:

  Matrix of events as produced by dd_sim:  
    
  - the first column is the time at which a species is born in Mya  
  - the second column is the label of the parent of the species;
  positive and negative values indicate whether the species belongs to
  the left or right crown lineage  
  - the third column is the label of the daughter species itself;
  positive and negative values indicate whether the species belongs to
  the left or right crown lineage  
  - the fourth column is the time of extinction of the species; if the
  fourth element equals -1, then the species is still extant.

- dropextinct:

  Sets whether the phylogeny should drop species that are extinct at the
  present

## Value

- phy :

  A phylogeny of the phylo type

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## Author

Rampal S. Etienne

## Examples

``` r
sim = dd_sim(c(0.2,0.1,20),10)
phy = L2phylo(sim$L)
plot(phy)

```
