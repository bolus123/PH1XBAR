\name{getCC}
\alias{getCC}
\title{get Phase I corrected charting constant}
 \usage{
getCC(
    m
    ,nu = m - 1
    ,FAP = 0.1
    ,off.diag = -1/(m - 1)
    ,var.est = 'MS'
    ,ub.option = TRUE
    ,maxiter = 10000
    ,ub.lower = 1e-6
    ,indirect.interval = c(1, 7)
    ,indirect.subdivisions = 100L
    ,indirect.tol = .Machine$double.eps^0.25
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
  \item{indirect.interval}{ the searching range of charting constants based on the indirect method }
  \item{indirect.subdivisions}{ the parameter of integration for the indirect method. }
  \item{indirect.tol}{ the tolerate level of convergence of integration for the indirect method. }

}
\description{
Obtain a corrected charting constant. If var.est is MS, the computation is carried out by the direct method. If var.est is MR, the computation is carried out by the indirect method. Notice that nu only works on when var.est is MS. When var.est is MR, nu will be recalculated.
}
\examples{
%\dontrun{
# MS-based estimator involved
getCC(m = 30, nu = 29, var.est = 'MS')

# MR-based estimator involved
getCC(m = 30, var.est = 'MR')
%}
}
