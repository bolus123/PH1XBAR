\name{PH1XBAR}
\alias{PH1XBAR}
\title{Build Phase I X-bar control chart with a corrected charting constant}
\usage{
PH1XBAR(
	X
	,cc = NULL
	,FAP0 = 0.1
	,var.est = c('S', 'MR')
	,ub.option = TRUE
	,method = c('exact', 'BA')
	,plot.option = TRUE
	,interval = c(1, 5)
	,nsim = 10000 
  	,verbose = FALSE
)
}
\arguments{
  \item{X}{ input and it must be a matrix }
  \item{cc}{ nominal Phase I charting constant. If this is given, the function will not recompute the charting constant. }
  \item{FAP0}{ nominal false Alarm Probabilty in Phase 1 }
  \item{var.est}{ 'S' - use mean-square-based estimator, 'MR' - use moving-range-based estimator}
  \item{ub.option}{ TRUE - the standard deviation estimator corrected by a unbiasing constant. For MS, it is c4 and for MR, it is d2. FALSE - no unbiasing constant }
  \item{method}{ 'exact' - calculate results using the exact method, 'BA' - calculate results using the Bonfferoni approximation }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{interval}{ searching range of charting constants for the exact method }
  \item{nsim}{ number of simulation for the exact method }
  \item{verbose}{ print diagnostic information about FAP0 and the charting constant during the simulations for the exact method }
}
\value{
Object of type list. A named list with the following objects:
  \item{CL}{Object type double - central line}
  \item{var.est}{Object type double - variance estimate}
  \item{ubCons}{Object type double - unbiasing constant}
  \item{cc}{Object type double - charting constant}
  \item{m}{Object type integer - number of observations}
  \item{nu}{Object type integer - degrees of freedom}
  \item{lambda}{Object type integer - chi-squared unbiasing constant}
  \item{LCL}{Object type double - lower charting limit}
  \item{UCL}{Object type double - upper charting limit}
  \item{CS}{Object type double - charting statistic}
}
\description{
Build a Phase I Shewhart control chart for the variance components model if the data are subgrouped or for the basic Shewhart model if the data are individual. The charting constant is correted by this approach.
}
\examples{

set.seed(12345)

# load the data in the package as an example
data(grinder_data)

# An example using a false alarm probability of 0.05, and 10 simulations
PH1XBAR(grinder_data, FAP0 = 0.05, nsim=10, verbose=TRUE)

}
