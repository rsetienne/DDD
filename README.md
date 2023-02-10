# DDD: Diversity-Dependent Diversification


[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/DDD)](https://cran.r-project.org/package=DDD)
[![](http://cranlogs.r-pkg.org/badges/grand-total/DDD)]( https://CRAN.R-project.org/package=DDD)
[![](http://cranlogs.r-pkg.org/badges/DDD)](https://CRAN.R-project.org/package=DDD)

Branch|[GitHub Actions](https://github.com/rsetienne/DDD/actions)|[Codecov](https://www.codecov.io)
---|---|---
`master`|[![Build Status](https://github.com/rsetienne/DDD/workflows/R-CMD-check/badge.svg?branch=master)](https://github.com/rsetienne/DDD/actions)|[![codecov.io](https://codecov.io/github/rsetienne/DDD/coverage.svg?branch=master)](https://codecov.io/github/rsetienne/DDD/branch/master)
`develop`|[![Build Status](https://github.com/rsetienne/DDD/workflows/R-CMD-check/badge.svg?branch=develop)](https://github.com/rsetienne/DDD/actions)|[![codecov.io](https://codecov.io/github/rsetienne/DDD/coverage.svg?branch=develop)](https://codecov.io/github/rsetienne/DDD/branch/develop)

Implements maximum likelihood and bootstrap methods based on
the diversity-dependent birth-death process to test whether
speciation or extinction are diversity-dependent, under various
models including various types of key innovations.

Also contains functions to simulate the diversity-dependent
process.

## Installing DDD

### From CRAN

From within R, do:

``` r
install.packages("DDD")
```

### From GitHub

Install `DDD` from this GitHub repository by running:

``` r
install.packages("remotes")
remotes::install_github("rsetienne/DDD")
```

## References
* Etienne et al. 2012, Proc. Roy. Soc. B 279: 1300-1309, https://doi.org/10.1098/rspb.2011.1439
* Etienne & Haegeman 2012, Am. Nat. 180: E75-E89, https://doi.org/10.1086/667574
* Etienne et al. 2016. Meth. Ecol. Evol. 7: 1092-1099, https://doi.org/10.1111/2041-210X.12565
* Laudanno et al. 2021. Syst. Biol. 70: 389–407, https://doi.org/10.1093/sysbio/syaa048
