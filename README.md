# Overview
The purpose of PH1XBAR is to build a Phase I Shewhart X-bar control chart with vairance components model in R.

# Installation

## Install from GitHub

PH1XBAR is still under development, so we recommend users installing it through Github as follows

``` r
install.packages("devtools")
devtools::install_github("bolus123/PH1XBAR")
```

Note that for Windows users,  Rtools may need to be installed in advance.  The detailed instruction is introduced: https://cran.r-project.org/bin/windows/Rtools/

For Mac and Linux users, please follow the instruction: https://www.r-project.org/nosvn/pandoc/devtools.html

## Install from local

Users can also download our release on Github and then install it from your local path as follows
``` r
install.packages('path_to_file', repos = NULL, type="source")
```


# Usage

Before using any functions, PH1XBAR may need to be loaded into R

``` r
library(PH1XBAR)
```

PH1XBAR provides a function to build Phase I X-bar chart with variance components model as follows

``` r
data(grinder_data)
PH1XBAR(grinder_data)
```

Notice that the variance estimator in the control chart must be S or MR. Also, PH1XBAR provides a function to get the corrected charting constant as follows

``` r
# S-based estimator involved
getCC(FAP0 = 0.1, m = 30, var.est = 'S')

# MR-based estimator involved
getCC(FAP0 = 0.1, m = 30, var.est = 'MR')
```

More details are on the manual.
