# Takes samples in the usual manner

The standard sample function in R samples from n numbers when x = n.
This is unwanted behavior when the size of the vector to sample from
changes dynamically. This is corrected in sample2

## Usage

``` r
sample2(x, size, replace = FALSE, prob = NULL)
```

## Arguments

- x:

  A vector of one or more elements

- size:

  A non-negative integer giving the number of items to choose.

- replace:

  Should sampling be with replacement?

- prob:

  A vector of probability weights for obtaining the elements of the
  vector being sampled.

## Value

- sam:

  A vector of length `size` that is sampled from `x`.

## Author

Rampal S. Etienne

## Examples

``` r
sample(x = 10,size = 5,replace = TRUE)
#> [1]  1 10  8  8  6
sample2(x = 10,size = 5,replace = TRUE)
#> [1] 10 10 10 10 10
```
