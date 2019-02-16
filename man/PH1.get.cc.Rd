\name{PH1.get.cc}
\alias{PH1.get.cc}
\title{PH1.get.cc}
 \usage{
PH1.get.cc(
    k
    , nu
    , FAP = 0.1
    , off.diag = -1/(k - 1)
    , c4.option = TRUE
    , maxiter = 10000
    , method = 'direct'
    , indirect.interval = c(1, 7)
    , indirect.subdivisions = 100L
    , indirect.tol = .Machine$double.eps^0.25
)
}
\arguments{
  \item{k}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{FAP}{ false Alarm Probabilty in Phase 1 }
  \item{off.diag}{ value of equal-correlation. Default setting is -1/(m - 1) }
  \item{c4.option}{ TRUE - the standard deviation estimator corrected by c4. FALSE - no c4 involved }
  %\item{alternative}{ '2-sided' - equal-limits applied in the calculation; '1-sided' - only upper limit involved }
  \item{maxiter}{ maximum iteration of computation in mvtnorm }
  \item{method}{ 'direct' - charting constants computed by multivariate t distribution. 'indirect' - charting constants computed by multivariate normal distribution }
  \item{indirect.interval}{ the searching range of charting constants based on the indirect method }
  \item{indirect.subdivisions}{ the parameter of integration for the indirect method. Please see __integrate__ }
  \item{indirect.tol}{ the tolerate level of convergence of integration for the indirect method. Please see __integrate__ }

}
\description{
Obtain a charting constant by the direct method or the indirect method.
}
\examples{
PH1.get.cc(
    k = 10
    ,nu = 40
)
}
