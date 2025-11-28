# Function to convert phylogeny to a table with speciation and extinction events

Converting a phylogeny to a table with speciation and extinction events

## Usage

``` r
phylo2L(phy)
```

## Arguments

- phy:

  A phylogeny of the phylo type

## Value

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

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## Author

Liang Xu

## Examples

``` r
sim = dd_sim(c(0.2,0.1,20),10)
phy = sim$tas
L = phylo2L(phy)
phy2 = L2phylo(L, dropextinct = FALSE)
graphics::par(mfrow = c(1,3))
graphics::plot(phy)
graphics::plot(phy2)
graphics::plot(L2phylo(sim$L, dropextinct = FALSE))

```
