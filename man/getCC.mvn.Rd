\name{getCC.mvn}
\alias{getCC.mvn}
\title{get Phase I corrected charting constant using the indirect method}
 \usage{
getCC.mvn(
    m
    ,nu = m - 1
    ,FAP = 0.1
    ,off.diag = -1/(m - 1)
    ,var.est = 'MS'
    ,ub.option = TRUE
    ,maxiter = 10000
    ,ub.lower = 1e-6
    ,interval = c(1, 7)
    ,subdivisions = 2000
    ,tol = 1e-2
)
}
\arguments{
  \item{m}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{FAP}{ false Alarm Probabilty in Phase 1 }
  \item{off.diag}{ value of equal-correlation. Default setting is -1/(m - 1) }
  \item{var.est}{ 'MS' - use mean-square-based estimator, 'MR' - use moving-range-based estimator}
  \item{ub.option}{ TRUE - the standard deviation estimator corrected by a unbiasing constant. For MS, it is c4 and for MR, it is d2. FALSE - no unbiasing constant }
  \item{maxiter}{ maximum iteration of computation in mvtnorm }
  \item{ub.lower}{ the lower bound of searching unbiasing constant for MR }
  \item{interval}{ the searching range of charting constants based on the indirect method }
  \item{subdivisions}{ the parameter of integration for the indirect method. }
  \item{tol}{ the tolerate level of convergence of integration for the indirect method. }

}
\description{
Obtain a corrected charting constant using the indirect method. Notice that this method is significantly slower than the direct method. We recommend users using the direct method instead if var.est is MS. Also, nu only works on when var.est is MS. When var.est is MR, nu will be recalculated.
}
\examples{
%\dontrun{
# MS-based estimator involved
getCC.mvn(m = 30, nu = 29, var.est = 'MS')

# MR-based estimator involved
getCC.mvn(m = 30, var.est = 'MR')
%}
}
