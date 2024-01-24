
<!-- README.md is generated from README.Rmd. Please edit that file -->

# validatEHR

validatEHR is a working version of the R package for the paper,
“Integrating EHR Subject to Outcome Misclassification with Validation
Data Subject to Selection Bias.” The package implements two approaches
for estimating the average treatment effect (ATE) when validation data
containing gold standard outcomes are available.

## Installation

You can install the development version of validatEHR from …

``` r
# install.packages("devtools")
devtools::install_github("jshen650/validatEHR")
```

## Example

This is a basic example which shows how to generate a data set from one
of the simulation scenarios described in the paper and apply one of the
two approaches described in the paper that leverages information from a
validation data set containing gold standard outcomes.

Data can be simulated with the `genDat` function. Data set size,
sensitivity of the outcomes, (1-specificity) of the outcomes, type of
validation data (e.g. “SRS” or “non-probability”), parameter for
validation data set size, and seed number can all be specified. The
corresponding arguments are “nA”, “p11”, “p10”, “valType”, “alphaInt”,
and “seedNum”.

``` r
## Simulate data following paper specifications - say n=5000, with non-probability validation sample of size nV~850

example_dat <- validatEHR::genDat(nA=5000, p11=0.67, p10=0.24, valType="non-probability", seedNum=215)
str(example_dat)
#> 'data.frame':    5000 obs. of  9 variables:
#>  $ y.A    : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ yStar  : int  1 0 0 0 0 0 1 0 0 1 ...
#>  $ B.index: int  0 0 0 0 0 0 1 0 0 0 ...
#>  $ t.A    : int  1 1 1 0 1 0 1 0 1 1 ...
#>  $ X1     : num  0.337 -1.045 1.471 1.509 -0.781 ...
#>  $ X2     : num  0.1402 -0.4261 -0.1201 0.0767 -1.0867 ...
#>  $ X3     : num  0.272 -1.644 -0.986 -0.159 0.123 ...
#>  $ X4     : num  1.4654 0.893 0.2559 -1.058 -0.0128 ...
#>  $ X5     : num  0.2721 0.4022 0.5815 -0.4815 0.0927 ...
```
