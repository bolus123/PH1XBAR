\name{getCC.ARMA}
\alias{getCC.ARMA}
\title{get Phase I corrected charting constant with an ARMA model}
 \usage{
getCC.ARMA(
	fap0 = 0.1 
  	,interval = c(1, 4)
  	,n = 50 
  	,order = c(1, 0, 0)
  	,phi.vec = 0.5
  	,theta.vec = NULL
  	,case = 'U'
  	,method = 'Method 3'
  	,nsim.coefs = 100
  	,nsim.process = 1000
  	,burn.in = 50
  	,sim.type = 'Matrix'
  	,logliktol = 1e-2
	,verbose = FALSE
)
}
\arguments{
  \item{fap0}{ nominal false Alarm Probabilty in Phase 1 }
  \item{interval}{ searching range of charting constants for the exact method }
  \item{n}{ number of observations }
  \item{order}{ order for ARMA model }
  \item{phi.vec}{ given vectors of autoregressive parameters for ARMA models }
  \item{theta.vec}{ given vectors of moving-average parameters for ARMA models }
  \item{case}{ known or unknown case.  When case = 'U', the parameters are estimated }
  \item{method}{ estimation method for the control chart. When method = 'Method 3' is maximum likehood estimations plus method of moments. Other options are 'Method 1' which is pure MLE and 'Method 2' which is pure CSS. }
  \item{nsim.coefs}{ number of simulation for coeficients.  It is functional when double.sim = TRUE. }
  \item{nsim.process}{ number of simulation for ARMA processes }
  \item{burn.in}{ number of burn-ins.  When burn.in = 0, the ECM gets involved.  When burn.in is large enough, the ACM gets involved. }
  \item{sim.type}{ type of simulation.  When sim.type = 'Matrix', the simulation is generated using matrix computation.  When sim.type = 'Recursive', the simulation is based on a recursion. }
  \item{logliktol}{ convergence tolerance for the log likelihood }
  \item{verbose}{ print diagnostic information about fap0 and the charting constant during the simulations for the exact method }
}
\value{
Object type double. The corrected charting constant. 
}
\description{
Obtain a corrected charting constant.  
}
\examples{

set.seed(12345)

# Calculate the charting constant using fap0 of 0.05, and 50 observations
getCC.ARMA(fap0=0.05, n=50, nsim.coefs=10, nsim.process=10)

}