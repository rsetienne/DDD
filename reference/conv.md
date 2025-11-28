# Function to do convolution of two vectors

Convolution of two vectors

## Usage

``` r
conv(x, y)
```

## Arguments

- x:

  first vector

- y:

  second vector

## Value

vector that is the convolution of x and y

## References

\- Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, doi:
10.1098/rspb.2011.1439  
- Etienne, R.S. & B. Haegeman 2012. Am. Nat. 180: E75-E89, doi:
10.1086/667574

## Author

Rampal S. Etienne

## Examples

``` r
conv(1:10,1:10)
#>  [1]   1   4  10  20  35  56  84 120 165 220 264 296 315 320 310 284 241 180 100
```
