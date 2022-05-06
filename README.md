[![R build
status](https://github.com/bolus123/PH1XBAR/workflows/R-CMD-check/badge.svg)](https://github.com/bolus123/PH1XBAR/actions)
[![](https://www.r-pkg.org/badges/version/PH1XBAR)](https://www.r-pkg.org/pkg/PH1XBAR)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/PH1XBAR)](https://www.r-pkg.org/pkg/PH1XBAR)

# Overview
The purpose of PH1XBAR is to build three types of Phase I Shewhart control charts:  
1. Phase I Shewhart X-bar control chart with a balanced one-way random effects model (doi:10.1002/qre.2793). 

2. Phase I Shewhart individual control chart for the iid case (doi:10.1080/08982112.2021.1878220). 

3. Phase I individual control chart with an ARMA model.  

# Installation

## Install from CRAN

PH1XBAR is published on CRAN, so we recommend users installing it in a regular way as follows

``` r
install.packages("PH1XBAR")
```

## Install from GitHub

PH1XBAR is still under development, so if users are more interested in the experimental version, there is an alternative installation through Github as follows

``` r
install.packages("devtools")
devtools::install_github("bolus123/PH1XBAR")
```

Note that for Windows users,  Rtools may need to be installed in advance.  Please choose the right version of Rtools which is corresponding to your R and Rstudio.  The detailed instruction is introduced: https://cran.r-project.org/bin/windows/Rtools/

For Mac and Linux users, please follow the instruction: https://www.r-project.org/nosvn/pandoc/devtools.html

## Install from local

Users can also download our release, PH1XBAR_x.y.z.tar.gz, from our homepage on CRAN or Github and then install it from your local path as follows
``` r
install.packages('path_to_file/PH1XBAR_x.y.z.tar.gz', repos = NULL, type="source")
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
getCC.XBAR(FAP0 = 0.1, m = 30, var.est = 'S')

# MR-based estimator involved
getCC.XBAR(FAP0 = 0.1, m = 30, var.est = 'MR')
```


PH1XBAR provides a function to build Phase I individual chart with an ARMA model as follows

``` r
data(preston_data)

# using the default setting whose FAP0 = 0.1
PH1ARMA(preston_data)

# using known parameters with FAP0 = 0.1
PH1ARMA(preston_data, case = 'K')
```

PH1XBAR provides a function to get the corrected charting constant for the ARMA model as follows

``` r
# Calculate the charting constant using FAP0 of 0.05, and 50 observations
getCC.ARMA(FAP0=0.05, n=50)
```

More details are on the manual.
