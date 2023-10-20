\name{PH1XBAR-package}
\alias{PH1XBAR-package}
\docType{package}
\title{
Phase I Shewhart X-bar Control Chart
}
\description{
The utility of this package is in building a Shewhart-type control chart based on new methods for subgrouped and individual data. The Phase I chart is based on the multivariate normal/t or ARMA process.
}
\details{
\tabular{ll}{
Package: \tab PH1XBAR\cr
Type: \tab Package\cr
% Version: \tab 0.11.1\cr
Date: \tab 2023-10-20\cr
License: \tab GPL (>= 3)\cr
}
Function 'PH1XBAR' builds a Phase I Shewhart X-bar/individual control chart with a correct charting constant.
Function 'getCC.XBAR' gets a charting constant for the data with a balanced one-way random effects model.
Function 'PH1ARMA' builds a Phase I individual control chart with an ARMA model using a correct charting constant.
Function 'getCC.ARMA' gets a charting constant for the data with an ARMA model.
}
\author{ Yuhui Yao, Subha Chakraborti, Tyler Thomas, Xin Yang, Jason Parton.

Maintainer: Yuhui Yao <yyao17@crimson.ua.edu>
}
\references{
Champ, C.W., and Jones, L.A. (2004) Designing Phase I X-bar charts with small sample sizes. Quality and Reliability Engineering International. 20(5), 497-510
Yao, Y., Hilton, C.W., and Chakraborti, S. (2017) Designing Phase I Shewhart X-bar charts: Extended tables and software. Quality and Reliability Engineering International. 33(8), 2667-2672.
Yao, Y., and Chakraborti, S. (2021). Phase I monitoring of individual normal data: Design and implementation. Quality Engineering, 33(3), 443-456.
Yao, Y., and Chakraborti, S. (2021). Phase I process monitoring: The case of the balanced one-way random effects model. Quality and Reliability Engineering International, 37(3), 1244-1265.
}
\examples{

# Build a Phase I basic Shewhart control chart
data(grinder_data)
PH1XBAR(grinder_data, nsim=10)

# Build a Phase I individual control chart with an ARMA model
data(preston_data)
PH1ARMA(preston_data, nsim.process=10, nsim.coefs=10)

}
