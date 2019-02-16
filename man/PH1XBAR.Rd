\name{PH1XBAR}
\alias{PH1XBAR}
\title{PH1XBAR}
\usage{
PH1XBAR(
	X
  ,c.i = NULL
	,FAP = 0.1
	,off.diag = -1/(k - 1)
	#,alternative = '2-sided'
  ,model = 'ANOVA-based'
  ,c4.option = TRUE
	,plot.option = TRUE
	,maxiter = 10000
	,method = 'direct'
	,indirect.interval = c(1, 7)
	,indirect.subdivisions = 100L
	,indirect.tol = .Machine$double.eps^0.25
)
}
\arguments{
  \item{X}{ input and it must be a matrix }
  \item{c.i}{ a nominal Phase I charting constant. If this is given, the function will not recompute the charting constant. }
  \item{FAP}{ false Alarm Probabilty in Phase 1 }
  \item{off.diag}{ value of equal-correlation. Default setting is -1/(m - 1) }
  %\item{alternative}{ '2-sided' - equal-limits applied in the calculation; '1-sided' - only upper limit involved }
  \item{model}{ 'ANOVA-based' - degrees of freedom equal to k - 1 and estimating the within-batch standard deviation. 'basic' - degree of freedo equal to k * (n - 1) and estimating the pooled standard deviation }
  \item{c4.option}{ TRUE - the standard deviation estimator corrected by c4. FALSE - no c4 involved }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{maxiter}{ maximum iteration of computation in mvtnorm }
  \item{method}{ 'direct' - charting constants computed by multivariate t distribution. 'indirect' - charting constants computed by multivariate normal distribution }
  \item{indirect.interval}{ the searching range of charting constants based on the indirect method }
  \item{indirect.subdivisions}{ the parameter of integration for the indirect method. Please see __integrate__ }
  \item{indirect.tol}{ the tolerate level of convergence of integration for the indirect method. Please see __integrate__ }
}
\description{
Build a Phase I Shewhart X-bar control chart by the basic model or the ANOVA-based model, and the direct method or the indirect method.
}
\examples{
X <- PH1XBAR.data()
PH1XBAR(X)
}
