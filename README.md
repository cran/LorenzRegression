
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LorenzRegression

<!-- badges: start -->

[![R-CMD-check](https://github.com/AlJacq/LorenzRegression/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/AlJacq/LorenzRegression/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The LorenzRegression package proposes a toolbox to estimate, produce
inference on and interpret Lorenz regressions. These regressions are
used to determine the explanatory power of a set of covariates on the
inequality of a response variable.

## Installation

You can install the released version of LorenzRegression from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("LorenzRegression")
```

You can install the development version of this package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AlJacq/LorenzRegression")
```

## What’s new

### LorenzRegression 2.0.0

- **Function Structure Overhaul**: `Lorenz.Reg` now acts as a wrapper
  for fitting functions `Lorenz.GA`, `Lorenz.FABS`, and
  `Lorenz.SCADFABS`, returning objects of class `"LR"` or `"PLR"` with
  designated methods.

- **Enhanced Bootstrap and CV**: `Lorenz.boot` performs bootstrap
  calculations and out-of-bag score computation, while `PLR.CV` handles
  cross-validation for tuning parameter selection.

- **Methods**: New methods `fitted`, `explainedIneq`, and `autoplot` for
  `"LR"` and `"PLR"` objects. Method availability is documented in the
  `Lorenz.Reg` help page with individual help pages for each method.
