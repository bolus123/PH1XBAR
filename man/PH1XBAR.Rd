\name{PH1XBAR}
\alias{PH1XBAR}
\title{Build Phase I X-bar control chart with a corrected charting constant}
\usage{
PH1XBAR(
    X
    ,cc = NULL
    ,FAP0 = 0.1
    ,off.diag = -1/(m - 1)
    ,var.est = c('S', 'MR')
    ,ub.option = TRUE
	,method = c('exact', 'BA')
    ,plot.option = TRUE
    ,interval = c(1, 5)
    ,nsim = 10000 
	,seed = 12345
)
}
\arguments{
  \item{X}{ input and it must be a matrix }
  \item{cc}{ nominal Phase I charting constant. If this is given, the function will not recompute the charting constant. }
  \item{FAP0}{ nominal false Alarm Probabilty in Phase 1 }
  \item{off.diag}{ value of equal-correlation. Default setting is -1/(m - 1) }
  \item{var.est}{ 'S' - use mean-square-based estimator, 'MR' - use moving-range-based estimator}
  \item{ub.option}{ TRUE - the standard deviation estimator corrected by a unbiasing constant. For MS, it is c4 and for MR, it is d2. FALSE - no unbiasing constant }
  \item{method}{ 'exact' - calculate results using the exact method, 'BA' - calculate results using the Bonfferoni approximation }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{interval}{ searching range of charting constants for the exact method }
  \item{nsim}{ number of simulation for the exact method }
  \item{seed}{ random seed for the exact method }
}
\description{
Build a Phase I Shewhart X-bar control chart with a variance component model. The charting constant is correted by this approach.
}
\examples{
%\dontrun{
data(grinder_data)
PH1XBAR(grinder_data, nsim = 1000)
%}
}
