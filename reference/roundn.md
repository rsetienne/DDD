# Rounds up in the usual manner

The standard round function in R rounds x.5 to the nearest even integer.
This is odd behavior that is corrected in roundn

## Usage

``` r
roundn(x, digits = 0)
```

## Arguments

- x:

  Number to be rounded

- digits:

  Sets the number of decimals in rounding.

## Value

- n:

  A number

## Author

Rampal S. Etienne

## Examples

``` r
round(2.5)
#> [1] 2
roundn(2.5)
#> [1] 3
round(3.5)
#> [1] 4
roundn(3.5)
#> [1] 4
round(2.65,digits = 1)
#> [1] 2.6
roundn(2.65,digits = 1)
#> [1] 2.7
round(2.75,digits = 1)
#> [1] 2.8
roundn(2.75,digits = 1)
#> [1] 2.8
```
