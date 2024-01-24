# validatEHR

validatEHR is a working version of the R package for the paper, "Integrating EHR Subject to Outcome Misclassification with Validation Data Subject to Selection Bias." The package implements two approaches for estimating the average treatment effect (ATE) when validation data containing gold standard outcomes are available.

## Installation

You can install the development version of validatEHR from ...

``` r
# install.packages("devtools")
devtools::install_github("jshen650/validatEHR")
```

## Example
This is a basic example which shows how to generate a data set from one of the simulation scenarios described in the paper and apply one of the two approaches described in the paper that leverages information from a validation data set containing gold standard outcomes.

Data can be simulated with the `genDat` function. Data set size, sensitivity of the outcomes, (1-specificity) of the outcomes, type of validation data (e.g. "SRS" or "non-probability"), parameter for validation data set size, and seed number can all be specified. The corresponding arguments are "nA", "p11", "p10", "valType", "alphaInt", and "seedNum".

``` r
## Simulate data following paper specifications - say n=5000, with non-probability validation sample of size nV~850

example_dat <- validatEHR::genDat(nA=5000, p11=0.67, p10=0.24, valType="non-probability", seedNum=215)
str(example_dat)

```



