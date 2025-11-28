# Carries out optimization (finding a minimum)

A wrapper to use several optimization routines, currently only 'simplex'
(a method adopted from Matlab, or 'subplex', from the R package
subplex). The function is called from several packages by the same
author.

## Usage

``` r
optimizer(
  optimmethod = "simplex",
  optimpars = c(1e-04, 1e-04, 1e-06, 1000),
  num_cycles = 1,
  fun,
  trparsopt,
  jitter = 0,
  ...
)
```

## Arguments

- optimmethod:

  The method to use for optimization, either 'simplex' or 'subplex'

- optimpars:

  Parameters of the optimization: 1) relative tolerance in function
  arguments, 2) relative tolerance in function value, 3) absolute
  tolerance in function arguments as well as the function value, 4)
  maximum number of iterations and 5) TRUE/FALSE flag to allow verbose
  output when using the simplex method, default is TRUE

- num_cycles:

  Number of cycles of the optimization. When set to Inf, the
  optimization will be repeated until the result is, within the
  tolerance, equal to the starting values, with a maximum of 10 cycles.

- fun:

  Function to be optimized

- trparsopt:

  Initial guess of the parameters to be optimized

- jitter:

  Perturbation of an initial parameter value when precisely equal to
  0.5; this is only relevant when subplex is chosen. The default value
  is 0, so no jitter is applied. A recommended value when using it is
  1E-5.

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
