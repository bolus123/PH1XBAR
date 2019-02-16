\name{PH1AND2XBAR-package}
\alias{00_PH1XBAR-package}
\docType{package}
\title{
Shewhart X-bar Control Chart in Phase I and in Phase II
}
\description{
The utility of this package is in building a Shewhart-type control chart based on new methods for the charting statistics, subgroup means. In Phase I, the chart is based on multivariate normal/t distribution. In Phase II, the chart is based on the uncondtional method or the exceedance probability criterion.
}
\details{
\tabular{ll}{
Package: \tab PH1XBAR\cr
Type: \tab Package\cr
% Version: \tab 1.0\cr
Date: \tab 2018-01-05\cr
License: \tab GPL (>= 3)\cr
}
Function 'PH1XBAR' builds a control chart by multivariate normal distribution.
Function 'get.cc' gets a charting constant by multivariate normal distribution.
}
\author{ Yuhui Yao, Subha Chakraborti.

Maintainer: Yuhui Yao <yyao17@crimson.ua.edu>
}
\references{
Champ CW, Jones LA. Designing Phase I X-bar charts with small sample sizes. Quality and Reliability Engineering International. 2004;20:497-510
Yao Y, Hilton CW, Chakraborti S. Designing Phase I Shewhart X-bar charts: Extended tables and software. Quality and Reliability Engineering International. 2017;33:2667–2672
Yao Y, Chakraborti S. Shewhart-Type Charts in Nonstandard Situations. 2018
}
\examples{
%\dontrun{
# Build a Phase I basic Shewhart control chart
X <- PH1XBAR.data()
PH1XBAR(X, model = 'basic')


# Build a Phase I ANOVA-based Shewhart control chart
Y <- PH2XBAR.data()
PH1XBAR(Y, model = 'ANOVA-based')

# Build a Phase II ANOVA-based Shewhart control chart
set.seed(1234)
Z<- matrix(rnorm(30 * 5, 245.1, sqrt(20.74) + sqrt(3.5)), nrow = 30, ncol = 5)
PH2XBAR(
    Z
    ,PH1.info = list(
    				X = Y
    				, mu = NULL
    				, sigma = NULL
    				, k = NULL
    				, n = NULL
    				, model = "ANOVA-based"
    			)
    ,c.ii.info = list(
    				c.ii = NULL
    				, method = 'both'
    				, ARL0 = 370
    				, p = 0.05
    				, eps = 0.1
    				, interval.c.ii.UC = c(1, 3.2)
    				, interval.c.ii.EPC = c(1, 10)
    				, UC.tol = .Machine$double.eps^0.25
            , EPC.tol = .Machine$double.eps^0.25
    			)
)

%}
}
