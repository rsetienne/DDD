# Transforming parameters from -Inf to Inf into parameters from -1 to 1

Function to transform pars in a way that is more useful for
optimization: trpars \<- sign(pars) \* pars/(sign(pars) + pars);

## Usage

``` r
transform_pars(pars)
```

## Arguments

- pars:

  Parameters to be transformed

## Value

Transformed parameters

## Author

Rampal S. Etienne
