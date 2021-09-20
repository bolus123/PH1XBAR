\name{PH1ARMA}
\alias{PH1ARMA}
\title{Build Phase I individual control chart with an ARMA model using a corrected charting constant}
\usage{
PH1ARMA(
	X
	,cc = NULL
	,FAP0 = 0.1
	,order = NULL
	,plot.option = TRUE
	,interval = c(1, 4)
	,case = 'U'
  	,method = 'Method 3'
	,nsimCoefs = 100
	,nsimProcess = 1000
	,burnIn = 50
	,simType = 'Matrix'
  	,logliktol = 1e-2
)
}
\arguments{
  \item{X}{ input and it must be a vector }
  \item{cc}{ nominal Phase I charting constant. If this is given, the function will not re-compute the charting constant. }
  \item{FAP0}{ nominal false Alarm Probabilty in Phase I }
  \item{order}{ order for ARMA model }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{interval}{ searching range of charting constants for the exact method }
  \item{case}{ known or unknown case.  When case = 'U', the parameters are estimated }
  \item{method}{ estimation method for the control chart. When method = 'Method 3' is maximum likehood estimations plus method of moments. Other options are 'Method 1' which is pure MLE and 'Method 2' which is pure CSS. }
  \item{nsimCoefs}{ number of simulation for coeficients.  It is functional when double.sim = TRUE. }
  \item{nsimProcess}{ number of simulation for ARMA processes }
  \item{burnIn}{ number of burn-ins.  When burnIn = 0, the ECM gets involved.  When burnIn is large enough, the ACM gets involved. }
  \item{simType}{ type of simulation.  When simType = 'Matrix', the simulation is generated using matrix computation.  When simType = 'Recursive', the simulation is based on a recursion. }
  \item{logliktol}{ convergence tolerance for the log likelihood }
}
\description{
Build a Phase I individual control chart with an ARMA model. The charting constant is correted by this approach.
}
\examples{

# load the data in the package as an example
data(preston_data)

# set number of simulations
nsimProcess <- 10
nsimCoefs <- 10

# An example using the default setting whose FAP0 = 0.1
PH1ARMA(preston_data, nsimProcess = nsimProcess, nsimCoefs = nsimCoefs)

# When users get an error message about the size of matrix,
# the function needs to use the alternative simulation type as follows
PH1ARMA(preston_data, FAP0 = 0.05, 
	nsimProcess = nsimProcess, nsimCoefs = nsimCoefs, simType = 'Recursive')

}
