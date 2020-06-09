\name{getCC}
\alias{getCC}
\title{get Phase I corrected charting constant}
 \usage{
getCC(
	FAP0
	,m 
	,var.est = c('S', 'MR') 
	,ubCons = 1 
	,method = c('exact', 'BA')
	,interval = c(1, 5)
	,nsim = 10000 
	,seed = 12345
	,nu = m - 1 
	,lambda = 1
)
}
\arguments{
  \item{FAP0}{ nominal false Alarm Probabilty in Phase 1 }
  \item{m}{ number of subgroups }
  \item{var.est}{ 'S' - use mean-square-based estimator, 'MR' - use moving-range-based estimator}
  \item{ubCons}{ unbiasing constant }
  \item{method}{ 'exact' - calculate results using the exact method, 'BA' - calculate results using the Bonfferoni approximation }
  \item{interval}{ searching range of charting constants for the exact method }
  \item{nsim}{ number of simulation for the exact method }
  \item{seed}{ random seed for the exact method }
  \item{nu}{ degrees of freedom for the Bonfferoni approximation}
  \item{lambda}{ constant for the Bonfferoni approximation }

}
\description{
Obtain a corrected charting constant.  
}
\examples{
%\dontrun{
# S-based estimator involved
getCC(FAP0 = 0.1, m = 20, var.est = 'S')

# MR-based estimator involved
getCC(FAP0 = 0.05, m = 20, var.est = 'MR')

# Calculate the charting constant using the BA method
getCC(FAP0 = 0.05, m = 20, method = 'BA')

%}
}
