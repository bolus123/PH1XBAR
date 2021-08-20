\name{PH1XBAR-package}
\alias{00_PH1XBAR-package}
\docType{package}
\title{
Phase I Shewhart X-bar Control Chart
}
\description{
The utility of this package is in building a Shewhart-type control chart based on new methods for the charting statistics, subgroup means. The Phase I chart is based on multivariate normal/t distribution.
}
\details{
\tabular{ll}{
Package: \tab PH1XBAR\cr
Type: \tab Package\cr
% Version: \tab 1.0\cr
Date: \tab 2018-01-05\cr
License: \tab GPL (>= 3)\cr
}
Function 'PH1XBAR' builds a Phase I Shewhart X-bar control chart with a correct charting constant.
Function 'getCC' gets a charting constant for the data with a balanced one-way random effects model.
}
\author{ Yuhui Yao, Subha Chakraborti.

Maintainer: Yuhui Yao <yyao17@crimson.ua.edu>
}
\references{
Champ CW, Jones LA. Designing Phase I X-bar charts with small sample sizes. Quality and Reliability Engineering International. 2004;20:497-510
Yao Y, Hilton CW, Chakraborti S. Designing Phase I Shewhart X-bar charts: Extended tables and software. Quality and Reliability Engineering International. 2017;33:2667–2672
}
\examples{
\dontrun{

# Build a Phase I basic Shewhart control chart
data(grinder_data)
PH1XBAR(grinder_data)


}
}
