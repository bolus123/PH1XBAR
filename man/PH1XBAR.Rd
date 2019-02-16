\name{PH1XBAR}
\alias{PH1XBAR}
\title{Build Phase I X-bar control chart with a corrected charting constant}
\usage{
PH1XBAR(
    X
    ,c.i = NULL
    ,FAP = 0.1
    ,off.diag = -1/(m - 1)
    ,var.est = 'MS'
    ,ub.option = TRUE
    ,plot.option = TRUE
    ,maxiter = 10000
    ,ub.lower = 1e-6
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
  \item{var.est}{ 'MS' - use mean-square-based estimator, 'MR' - use moving-range-based estimator}
  \item{ub.option}{ TRUE - the standard deviation estimator corrected by a unbiasing constant. For MS, it is c4 and for MR, it is d2. FALSE - no unbiasing constant }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{maxiter}{ maximum iteration of computation in mvtnorm }
  \item{ub.lower}{ the lower bound of searching unbiasing constant for MR }
  \item{indirect.interval}{ the searching range of charting constants based on the indirect method }
  \item{indirect.subdivisions}{ the parameter of integration for the indirect method. }
  \item{indirect.tol}{ the tolerate level of convergence of integration for the indirect method. }
}
\description{
Build a Phase I Shewhart X-bar control chart with a variance component model. The charting constant is correted by this approach.
}
\examples{
%\dontrun{
data(grinder_data)
PH1XBAR(grinder_data)
%}
}
