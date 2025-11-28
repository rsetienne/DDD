# Sampling in which zero probabilities are removed

Sampling in which zero probabilities are removed

## Usage

``` r
rng_respecting_sample(x, size, replace, prob)
```

## Arguments

- x:

  either a vector of one or more elements from which to choose, or a
  positive integer. See ‘Details.’

- size:

  a non-negative integer giving the number of items to choose.

- replace:

  should sampling be with replacement?

- prob:

  a vector of probability weights for obtaining the elements of the
  vector being sampled.

## Value

a vector of length size with elements drawn from either x or from the
integers 1:x.

## Note

thanks to Pedro Neves for finding this feature in base::sample

## See also

See [`sample`](https://rdrr.io/r/base/sample.html) for more details

## Author

Richel J.C. Bilderbeek

## Examples

``` r
  # Number of draws
  n <- 1000
  
  # Do normal sampling
  set.seed(42)
  draws_1 <- DDD:::rng_respecting_sample(
    1:3, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0)
  )
  
  # Do a sampling with one element of probability zero
  set.seed(42)
  draws_2 <- DDD:::rng_respecting_sample(
    1:4, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0, 0.0)
  )
  testit::assert(sum(draws_2 == 4) == 0)
  testit::assert(draws_1 == draws_2)
  
  # Use base sampling will give different results,
  # as it results in different RNG values
  set.seed(42)
  draws_3 <- sample(
    1:4, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0, 0.0)
  )
  testit::assert(sum(draws_3 == 4) == 0)
  testit::assert(!all(draws_1 == draws_3))
  
```
