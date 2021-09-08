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
  \item{nu}{ degrees of freedom for the Bonfferoni approximation}
  \item{lambda}{ constant for the Bonfferoni approximation }
}
\value{
Object type double. The corrected charting constant. 
}
\description{
Obtain a corrected charting constant.  
}
\examples{

set.seed(12345)

# Calculate the charting constant using 10 simulations and mean-square-based estimator
getCC(FAP0=0.05, m=20, nsim=10, var.est='S')

# Calculate the charting constant using 10 simulations and moving-range-based estimator
getCC(FAP0=0.05, m=20, nsim=10, var.est='MR')

}
