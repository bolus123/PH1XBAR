# Overview
The purpose of PH1XBAR is to build a Phase I Shewhart X-bar control chart with vairance components model in R.

# Installation

## Install from GitHub

PH1XBAR is still under development, so we recommend users installing it through Github as follows

``` r
install.packages("devtools")
devtools::install_github("bolus123/PH1XBAR")
```

## Install from local

Users can also download our release on Github and then install it from your local path as follows
``` r
install.packages('path_to_file', repos = NULL, type="source")
```


# Usage

PH1XBAR provides a function to build Phase I X-bar chart with variance components model as follows

``` r
data(grinder_data)
PH1XBAR(grinder_data)
```

Notice that the variance estimator in the control chart must be MS or MR. Also, PH1XBAR provides a function to get the corrected charting constant as follows

``` r
# MS-based estimator involved
getCC(m = 30, nu = 29, var.est = 'MS')

# MR-based estimator involved
getCC(m = 30, var.est = 'MR')
```

More details are on the manual.
