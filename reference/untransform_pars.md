# Untransforming parameters from -1 to 1 into parameters from -Inf to Inf.

Function to untransform pars after optimization: pars \<- sign(trpars)
\* trpars/(sign(trpars) - trpars);

## Usage

``` r
untransform_pars(trpars)
```

## Arguments

- trpars:

  Parameters to be untransformed

## Value

Untransformed parameters

## Author

Rampal S. Etienne
