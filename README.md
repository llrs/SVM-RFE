
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mRFE

<!-- badges: start -->

[![R-CMD-check](https://github.com/llrs/SVM-RFE/workflows/R-CMD-check/badge.svg)](https://github.com/llrs/SVM-RFE/actions)
<!-- badges: end -->

The goal of mRFE is to provide a package for mSVM-RFE: (multiple)
Support Vector Machine Recursive Feature Elimination

This repository contains an [R](http://www.r-project.org) implementation
of the mSVM-RFE algorithm ([Duan et al.,
2005](http://www.ncbi.nlm.nih.gov/pubmed/16220686)), including the
option to cut the features by half each round (instead of one-by-one) if
there are many features.

Also included are tools for wrapping the feature ranking/selection
process in an external layer of cross-validation for obtaining unbiased
estimates of generalization error/accuracy (See [Ambroise et al.,
2002](http://www.ncbi.nlm.nih.gov/pubmed/11983868)).

## Installation

You can install the development version of mRFE from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("llrs/SVM-RFE")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(mRFE)
#> Loading required package: e1071
## basic example code
```

## References

Please use `citation(mRFE)` to cite the original work:

### SVM-RFE

An iterative algorithm that works backward from an initial set of
features. At each round it 1. fits a simple linear SVM, 2. ranks the
features based on their weights in the SVM solution, and 3. eliminates
the feature with the lowest weight.

### *Multiple* SVM-RFE

Extends this idea by using resampling techniques at each iteration to
stabilize the feature rankings. Here we use cross validation. The
mSVM-RFE paper is:
