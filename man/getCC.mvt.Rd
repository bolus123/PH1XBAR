\name{getCC.mvt}
\alias{getCC.mvt}
\title{get Phase I corrected charting constant using the direct method}
 \usage{
getCC.mvt(
    m
    ,nu = m - 1
    ,FAP = 0.1
    ,off.diag = -1/(m - 1)
    ,ub.option = TRUE
    ,maxiter = 10000
)
}
\arguments{
  \item{m}{ number of subgroups }
  \item{nu}{ degrees of freedom }
  \item{FAP}{ false Alarm Probabilty in Phase 1 }
  \item{off.diag}{ value of equal-correlation. Default setting is -1/(m - 1) }
  \item{ub.option}{ TRUE - the standard deviation estimator corrected by c4. FALSE - no unbiasing constant }
  \item{maxiter}{ maximum iteration of computation in mvtnorm. }
}
\description{
Obtain a corrected charting constant using the direct method.
}
\examples{
%\dontrun{
getCC.mvt(m = 10, nu = 9)
%}
}
