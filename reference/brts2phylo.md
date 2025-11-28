# Function to convert a set of branching times into a phylogeny with random topology This code is taken from the package TESS by Sebastian Hoehna, where the function is called tess.create.phylo

Converting a set of branching times to a phylogeny

## Usage

``` r
brts2phylo(times, root = FALSE, tip.label = NULL)
```

## Arguments

- times:

  Set of branching times

- root:

  When root is FALSE, the largest branching time will be assumed to be
  the crown age. When root is TRUE, it will be the stem age.

- tip.label:

  Tip labels. If set to NULL, the labels will be t1, t2, etc.

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
