# varEqChol
## Introduction

The `varEqChol` is a package that learns ordering from the Cholesky factor of inverse covariance matrix under assumption that vriances of error terms are equal. 

This document serves as an introduction of using the package.

The main function is `eqvarDAG_ch`, which takes a data matrix of the observations and returns the ordering of the variables. 

## Installation

To install the latest version from Github, use

```s
library(devtools)
devtools::install_github("adallak/varEqChol")
```
