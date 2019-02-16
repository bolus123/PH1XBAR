\name{getCC.mvn.MC}
\alias{getCC.mvn.MC}
\title{get Phase I corrected charting constant using the indirect method with Monte Carlo simulation}
 \usage{
getCC.mvn.MC(
    m
    ,nu = m - 1
    ,FAP = 0.1
    ,off.diag = -1/(m - 1)
    ,var.est = 'MS'
    ,ub.option = TRUE
    ,sim.X = 10000
    ,sim.Y = 10000
    ,maxiter = 10000
    ,ub.lower = 1e-6
    ,interval = c(1, 7)
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
  \item{sim.X}{ number of simulation of random samples from multivariate normal distribution}
  \item{sim.Y}{ number of simulation of random samples from chi-square distribution with nu degrees of freedom }
  \item{maxiter}{ maximum iteration of computation in mvtnorm }
  \item{ub.lower}{ the lower bound of searching unbiasing constant for MR }
  \item{interval}{ the searching range of charting constants based on the indirect method }
  \item{tol}{ the tolerate level of convergence of integration for the indirect method.}

}
\description{
Obtain a corrected charting constant using the indirect method with Monte Carlo simulation. Notice that this method is significantly slower and much more inaccurate than other methods. We recommend users using other methods except for special demands.
}
\examples{
%\dontrun{
# MS-based estimator involved
getCC.mvn.MC(m = 30, nu = 29, var.est = 'MS')

# MR-based estimator involved
getCC.mvn.MC(m = 30, var.est = 'MR')
%}
}
