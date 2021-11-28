\name{getcc.arma}
\alias{getcc.arma}
\title{get Phase I corrected charting constant with an ARMA model}
 \usage{
getcc.arma(
	FAP0 = 0.1 
  	,interval = c(1, 4)
  	,n = 50 
  	,order = c(1, 0, 0)
  	,phiVec = 0.5
  	,thetaVec = NULL
  	,case = 'U'
  	,method = 'Method 3'
  	,nsimCoefs = 100
  	,nsimProcess = 1000
  	,burnIn = 50
  	,simType = 'Matrix'
  	,logliktol = 1e-2
	,verbose = FALSE
)
}
\arguments{
  \item{FAP0}{ nominal false Alarm Probabilty in Phase 1 }
  \item{interval}{ searching range of charting constants for the exact method }
  \item{n}{ number of observations }
  \item{order}{ order for ARMA model }
  \item{phiVec}{ given vectors of autoregressive parameters for ARMA models }
  \item{thetaVec}{ given vectors of moving-average parameters for ARMA models }
  \item{case}{ known or unknown case.  When case = 'U', the parameters are estimated }
  \item{method}{ estimation method for the control chart. When method = 'Method 3' is maximum likehood estimations plus method of moments. Other options are 'Method 1' which is pure MLE and 'Method 2' which is pure CSS. }
  \item{nsimCoefs}{ number of simulation for coeficients.  It is functional when double.sim = TRUE. }
  \item{nsimProcess}{ number of simulation for ARMA processes }
  \item{burnIn}{ number of burn-ins.  When burnIn = 0, the ECM gets involved.  When burnIn is large enough, the ACM gets involved. }
  \item{simType}{ type of simulation.  When simType = 'Matrix', the simulation is generated using matrix computation.  When simType = 'Recursive', the simulation is based on a recursion. }
  \item{logliktol}{ convergence tolerance for the log likelihood }
  \item{verbose}{ print diagnostic information about FAP0 and the charting constant during the simulations for the exact method }
}
\value{
Object type double. The corrected charting constant. 
}
\description{
Obtain a corrected charting constant.  
}
\examples{

set.seed(12345)

# Calculate the charting constant using FAP0 of 0.05, and 50 observations
getcc.arma(FAP0=0.05, n=50, nsimCoefs=10, nsimProcess=10)

}
