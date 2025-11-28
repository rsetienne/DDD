# Carries out optimization using a simplex algorithm (finding a minimum)

Function to optimize target function using a simplex method adopted from
Matlab

## Usage

``` r
simplex(fun, trparsopt, optimpars, ...)
```

## Arguments

- fun:

  Function to be optimized

- trparsopt:

  Initial guess of the parameters to be optimized

- optimpars:

  Parameters of the optimization: 1) relative tolerance in function
  arguments, 2) relative tolerance in function value, 3) absolute
  tolerance in function arguments as well as the function value, 4)
  maximum number of iterations and 5) TRUE/FALSE flag to allow verbose
  output, default is TRUE

- ...:

  Any other arguments of the function to be optimimized, or settings of
  the optimization routine

## Value

- out:

  A list containing optimal function arguments (`par`, the optimal
  function value (`fvalues`) and whether the optimization converged
  (`conv`)

.

## Author

Rampal S. Etienne

## Examples

``` r
cat("No examples")
#> No examples
```
