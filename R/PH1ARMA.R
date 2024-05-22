#' Phase I individual control chart with an ARMA model
#' 
#' Build a Phase I individual control chart for the ARMA models. The charting constant is corrected by this approach.
#' @param X input and it must be a vector 
#' @param cc nominal Phase I charting constant. If this is given, the function will not re-compute the charting constant. 
#' @param fap0 nominal false Alarm Probabilty in Phase I 
#' @param order order for ARMA model 
#' @param plot.option - draw a plot for the process; FALSE - Not draw a plot for the process
#' @param interval searching range of charting constants for the exact method 
#' @param case known or unknown case.  When case = 'U', the parameters are estimated 
#' @param phi.vec vector of autoregressive coefficient(s).  When case = 'K', the vector needs to be provided with the length same as the first value in the order.  If autoregressive coefficents does not present, phi needs to be NULL
#' @param theta.vec vector of moving-average coefficient(s).  When case = 'K', the vector needs to be provided with the length same as the third value in the order.  If moving-average coefficents does not present, theta needs to be NULL
#' @param mu0 value of the IC process mean.  When case = 'K', the value needs to be provided.
#' @param sigma0 value of the IC process standard deviation.  When case = 'K', the value needs to be provided.
#' @param method estimation method for the control chart. When method = 'MLE+MOM' is maximum likehood estimations plus method of moments. Other options are 'MLE' which is pure MLE and 'CSS' which is pure CSS. 
#' @param nsim.coefs number of simulation for coeficients.
#' @param nsim.process number of simulation for ARMA processes 
#' @param burn.in number of burn-ins.  When burn.in = 0, the simulated process is assumed to be in the initial stage.  When burn.in is large enough, the simulated process is assumed to be in the stable stage. 
#' @param sim.type type of simulation.  When sim.type = 'Matrix', the simulation is generated using matrix computation.  When sim.type = 'Recursive', the simulation is based on a recursion. 
#' @param standardize Output standardized data instead of raw data 
#' @param verbose print diagnostic information about fap0 and the charting constant during the simulations for the exact method 
#' @return CL Object type double - central line
#' @return gamma Object type double - process variance estimate
#' @return cc Object type double - charting constant
#' @return order Object type integer - order for ARMA model
#' @return phi.vec Object type integer - values of autoregressors
#' @return theta.vec Object type integer - values of moving averages
#' @return LCL Object type double - lower charting limit
#' @return UCL Object type double - upper charting limit
#' @return CS Object type double - charting statistic
#' @references Yao, Y., Chakraborti, S., Yang, X., Parton, J., Lewis Jr, D., and Hudnall, M. (2023). Phase I control chart for individual autocorrelated data: application to prescription opioid monitoring. Journal of Quality Technology, 55(3), 302-317.
#' 
#' 
#' @export
#' @examples
#' # load the data in the package as an example
#' data(preston_data)
#' 
#' # set number of simulations
#' nsim.process <- 10
#' nsim.coefs <- 10
#' 
#' # An example using the default setting whose fap0 = 0.1
#' PH1ARMA(preston_data, nsim.process = nsim.process, nsim.coefs = nsim.coefs)
#' 
#' # When users get an error message about the size of matrix,
#' # the function needs to use the alternative simulation type as follows
#' PH1ARMA(preston_data, fap0 = 0.05, 
#' 	nsim.process = nsim.process, nsim.coefs = nsim.coefs, sim.type = 'Recursive')
#' 
PH1ARMA <- function(X, cc = NULL, fap0 = 0.05, order = c(1, 0, 0), plot.option = TRUE, interval = c(1, 4),
                    case = 'U', phi.vec = NULL, theta.vec = NULL, mu0 = NULL, sigma0 = NULL, method = 'MLE+MOM', nsim.coefs = 100, nsim.process = 1000, burn.in = 50, 
                    sim.type = 'Matrix', standardize=TRUE, verbose = FALSE) {

  if (!is.vector(X)) {
	if (dim(X)[1] == 1 | dim(X)[2] == 1) {
		X <- as.vector(X)
	} else {
		stop('X is not a vector, or a m x 1 or 1 x m matrix.')
	}
  } 

  n <- length(X)

  if (is.null(cc)) {

    if (case == 'U') {
      if ((method == "MLE")||(method == "MLE+MOM")) {
        model <- arima(X, order = order, method = "CSS-ML")
      } else {
        model <- arima(X, order = order, method = "CSS")
      }
      
      
      if (length(model$model$phi) > 0) {
        phi.vec <- model$model$phi
      } else {
        phi.vec <- NULL
      }
  
      if (length(model$model$theta) > 0) {
        theta.vec <- model$model$theta
      } else {
        theta.vec <- NULL
      }
  
      cc <- getCC.ARMA(fap0 = fap0, interval = interval, n, order = order, phi.vec = phi.vec, theta.vec = theta.vec, case = case,
        method = method, nsim.coefs = nsim.coefs, nsim.process = nsim.process, burn.in = burn.in, sim.type = sim.type, verbose = verbose)

    } else if (case == 'K') {
      cc <- getCC.ARMA(fap0 = fap0, interval = interval, n, order = order, phi.vec = phi.vec, theta.vec = theta.vec, case = case,
        method = method, nsim.coefs = nsim.coefs, nsim.process = nsim.process, burn.in = burn.in, sim.type = sim.type, verbose = verbose)

    }
    
    
  }

  if (order[2] > 0) {

    X <- diff(X, differences = order[2])

  }

  if (case == "U") {
    if (standardize){
      mu <- mean(X)
      gamma <- sd(X)
  
      stdX <- (X - mu) / gamma
  
      LCL <- -cc
      UCL <- cc
    }else{
      mu <- mean(X)
      gamma <- sd(X)
  
      stdX <- X
  
      LCL <- -cc * gamma + mu
      UCL <- cc * gamma + mu
    }
  } else if (case == "K") {
    
    if (standardize){
      mu <- mu0
      gamma <- sigma0
  
      stdX <- (X - mu) / gamma
  
      LCL <- -cc
      UCL <- cc
    }else{
      mu <- mu0
      gamma <- sigma0
  
      stdX <- X
  
      LCL <- -cc * gamma + mu
      UCL <- cc * gamma + mu
    }
    
  }
  
  

  if (plot.option == TRUE) {

    main.text <- paste('Phase I Individual Chart for fap0 =', fap0, 'with an ARMA model')

    plot(c(1, n), c(min(LCL, stdX), max(UCL, stdX)), xaxt = "n", xlab = 'Observation', ylab = 'Charting Statistic', type = 'n', main = main.text)

    axis(side = 1, at = 1:n)

    points((1 + order[2]):n, stdX, type = 'o')
    points(c(-1, n + 2), c(LCL, LCL), type = 'l', col = 'red')
    points(c(-1, n + 2), c(UCL, UCL), type = 'l', col = 'red')
    points(c(-1, n + 2), c(mu, mu), type = 'l', col = 'blue')
    text(round(n * 0.8), UCL, paste('UCL = ', round(UCL, 4)), pos = 1)
    text(round(n * 0.8), LCL, paste('LCL = ', round(LCL, 4)), pos = 3)

  }

  res <- list(CL = mu, gamma = gamma, cc = cc, order = order, phi.vec = phi.vec, theta.vec = theta.vec, LCL = LCL, UCL = UCL, CS = stdX)

  return(invisible(res))
}
