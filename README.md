
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

Data can be simulated with the `genDat` function. This function requires
specification of the following:  
\* n = data set size  
\* p11 = sensitivity of the outcomes  
\* p10 = (1-specificity) of the outcomes  
\* valType = type of validation data (e.g. “SRS” or “non-probability”)  
\* alphaInt = parameter for validation data set size; if not specified,
the default will generate size of \~850  
\* seedNum = seed number

``` r
## Simulate data following paper specifications - say n=5000, with non-probability validation sample of size nV~850

example_dat <- validatEHR::genDat(n=5000, p11=0.67, p10=0.24, valType="non-probability", seedNum=215)
str(example_dat)
#> 'data.frame':    5000 obs. of  9 variables:
#>  $ y    : int  0 0 0 0 0 0 0 0 0 0 ...
#>  $ yStar: int  1 0 0 0 0 0 1 0 0 1 ...
#>  $ inVal: int  0 0 0 0 0 0 1 0 0 0 ...
#>  $ trt  : int  1 1 1 0 1 0 1 0 1 1 ...
#>  $ X1   : num  0.337 -1.045 1.471 1.509 -0.781 ...
#>  $ X2   : num  0.1402 -0.4261 -0.1201 0.0767 -1.0867 ...
#>  $ X3   : num  0.272 -1.644 -0.986 -0.159 0.123 ...
#>  $ X4   : num  1.4654 0.893 0.2559 -1.058 -0.0128 ...
#>  $ X5   : num  0.2721 0.4022 0.5815 -0.4815 0.0927 ...
```

The estimate of the ATE and its variance using gold standard outcomes
from the validation data and silver standard outcomes from the
non-validation individuals can be obtained using `run_method2`, which
follows from applying Equation 2 in the paper. Alternatively, for using
all silver standard outcomes rather than just a subset, estimation of
the ATE and its variance can come from implementing `run_method4`. The
user also has the option to use the optimal choice of weight, $w$, which
leads to minimum variance for obtaining the final estimate of the ATE.

``` r
## estimate ATE and variance when using gold standard outcomes from validation sample
## and all silver standard outcomes

## specify relevant variables for treatment propensity model
varTrtMod = c("X1", "X2", "X3", "X4", "X5")

## specify relevant variables for validation sample selection model - 
## here, includes the same covariates as varTrtMod
varSelectMod = varTrtMod

example_ATE <- validatEHR::run_est4(example_dat, inVal=example_dat$inVal, varTrt="trt",
                                       varGold="y", varSilver="yStar", varTrtMod=varTrtMod, 
                                       varSelectMod = varSelectMod, opt="TRUE")

## display optimal choice of weight w, final estimate of ATE, final variance estimate of ATE,
# estimate of ATE using validation data, variance estimate of ATE using validation data,
# estimate of ATE using silver standard outcomes,
# and variance estimate of ATE using silver standard outcomes
example_ATE
#>          w    estATE       varATE estATE_gold  varATE_gold estATE_silver
#> 1 0.605865 0.1112225 0.0005467645   0.1067438 0.0008637394     0.1181072
#>   varATE_silver
#> 1   0.001295773
```
