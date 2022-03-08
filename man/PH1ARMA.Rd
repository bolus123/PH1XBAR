\name{PH1ARMA}
\alias{PH1ARMA}
\title{Build Phase I individual control chart with an ARMA model using a corrected charting constant}
\usage{
PH1ARMA(
	X
	,cc = NULL
	,fap0 = 0.1
	,order = NULL
	,plot.option = TRUE
	,interval = c(1, 4)
	,case = 'U'
	,method = 'Method 3'
	,nsim.coefs = 100
	,nsim.process = 1000
	,burn.in = 50
	,sim.type = 'Matrix'
	,logliktol = 1e-2
	,verbose = FALSE
	,max.p=1
	,max.q=0
	,max.d=0
  ,standardize=TRUE
)
}
\arguments{
  \item{X}{ input and it must be a vector }
  \item{cc}{ nominal Phase I charting constant. If this is given, the function will not re-compute the charting constant. }
  \item{fap0}{ nominal false Alarm Probabilty in Phase I }
  \item{order}{ order for ARMA model }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{interval}{ searching range of charting constants for the exact method }
  \item{case}{ known or unknown case.  When case = 'U', the parameters are estimated }
  \item{method}{ estimation method for the control chart. When method = 'Method 3' is maximum likehood estimations plus method of moments. Other options are 'Method 1' which is pure MLE and 'Method 2' which is pure CSS. }
  \item{nsim.coefs}{ number of simulation for coeficients.}
  \item{nsim.process}{ number of simulation for ARMA processes }
  \item{burn.in}{ number of burn-ins.  When burn.in = 0, the ECM gets involved.  When burn.in is large enough, the ACM gets involved. }
  \item{sim.type}{ type of simulation.  When sim.type = 'Matrix', the simulation is generated using matrix computation.  When sim.type = 'Recursive', the simulation is based on a recursion. }
  \item{logliktol}{ convergence tolerance for the log likelihood }
  \item{verbose}{ print diagnostic information about fap0 and the charting constant during the simulations for the exact method }
  \item{max.p}{ max.p for auto.arima }
  \item{max.q}{ max.q for auto.arima }
  \item{max.d}{ max.d for auto.arima }
  \item{standardize}{ Output z-scores instead of raw data }
}
\value{
Object of type list. A named list with the following objects:
  \item{CL}{Object type double - central line}
  \item{gamma}{Object type double - process variance estimate}
  \item{cc}{Object type double - charting constant}
  \item{order}{Object type integer - order for ARMA model}
  \item{phi.vec}{Object type integer - values of autoregressors}
  \item{theta.vec}{Object type integer - values of moving averages}
  \item{LCL}{Object type double - lower charting limit}
  \item{UCL}{Object type double - upper charting limit}
  \item{CS}{Object type double - charting statistic}
}
\description{
Build a Phase I individual control chart for the ARMA models. The charting constant is corrected by this approach.
}
\examples{

# load the data in the package as an example
data(preston_data)

# set number of simulations
nsim.process <- 10
nsim.coefs <- 10

# An example using the default setting whose fap0 = 0.1
PH1ARMA(preston_data, nsim.process = nsim.process, nsim.coefs = nsim.coefs)

# When users get an error message about the size of matrix,
# the function needs to use the alternative simulation type as follows
PH1ARMA(preston_data, fap0 = 0.05, 
	nsim.process = nsim.process, nsim.coefs = nsim.coefs, sim.type = 'Recursive')

}
