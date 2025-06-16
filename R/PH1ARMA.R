#' Phase I individual control chart with an ARMA model
#' 
#' Build a Phase I individual control chart for the ARMA models. The charting constant is corrected by this approach.
#' @param X input and it must be a vector (m by 1)
#' @param cc nominal Phase I charting constant. If this is given, the function will not re-compute the charting constant. 
#' @param fap0 nominal false Alarm Probabilty in Phase I 
#' @param order order for ARMA(p, q) model 
#' @param plot.option - draw a plot for the process; TRUE - Draw a plot for the process, FALSE - Not draw a plot for the process
#' @param interval searching range of charting constants for the exact method 
#' @param case known or unknown case.  When case = 'U', the parameters are estimated, when case = 'K', the parameters need to be input
#' @param phi.vec a vector of length p containing autoregressive coefficient(s).  When case = 'K', the vector must have a length equal to the first value in the order.  If no autoregressive coefficent presents, set phi.vec = NULL
#' @param theta.vec a vector of length q containing moving-average coefficient(s).  When case = 'K', the vector must have a length equal to the first value in the order.  If no moving-average coefficent presents, set theta.vec = NULL
#' @param mu0 value of the IC process mean.  When case = 'K', the value needs to be provided.
#' @param sigma0 value of the IC process standard deviation.  When case = 'K', the value needs to be provided.
#' @param method estimation method for the control chart. When method = 'MLE+MOM' is maximum likehood estimations plus method of moments. Other options are 'MLE' which is pure MLE and 'CSS' which is pure CSS. 
#' @param nsim.coefs number of simulation for coefficients.
#' @param nsim.process number of simulation for ARMA processes 
#' @param burn.in number of burn-ins.  When burn.in = 0, the simulated process is assumed to be in the initial stage.  When burn.in is sufficiently large (e.g., the default value of 50), the simulated process is assumed to have reached a stable state. 
#' @param sim.type type of simulation. When sim.type = 'Recursive', the simulation is generated recursively, as in the ARMA model. When sim.type = 'Matrix', the simulation is generated using the covariance matrix among observations, derived from the relationship between the ARMA coefficient(s) and the partial autocorrelation(s). Note that sim.type = 'Matrix' is primarily used as a proof of concept and is not recommended for practical use due to its high computational cost.    
#' @param transform type of transformation. When transform = 'none', no transformation is performed. When transform = 'boxcox', the Box-Cox transformation is used. When transform = 'yeojohnson', the Yeo-Johnson transformation is used. 
#' @param lambda parameter used in the Box-Cox or Yeo-Johnson transformation.
#' @param standardize Output standardized charting statistics instead of raw ones. When standardize = TRUE, the standardization is used.  When standardize = FALSE, the standardization is not performed.
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
PH1ARMA <- function(X, cc = NULL, fap0 = 0.05, order = c(1, 0), plot.option = TRUE, interval = c(1, 4),
                    case = 'U', phi.vec = NULL, theta.vec = NULL, mu0 = NULL, sigma0 = NULL, method = 'MLE+MOM', nsim.coefs = 100, nsim.process = 1000, burn.in = 50, 
                    sim.type = 'Recursive', transform = "none", lambda = 1, standardize=FALSE, verbose = FALSE) {
  
  neworder <- c(order[1], 0, order[2])
  
  if (!is.vector(X)) {
  	if (dim(X)[1] == 1 | dim(X)[2] == 1) {
  		X <- as.vector(X)
  	} else {
  		stop('X is not a vector, or a m x 1 or 1 x m matrix.')
  	}
  } 

  lambda1 <- NA
  if (transform == "boxcox") {
    lambda1 <- lambda
    X <- forecast::BoxCox(X, lambda1)
  } else if (transform == "yeojohnson") {
    lambda1 <- lambda
    X <- VGAM::yeo.johnson(X, lambda1)
  }
  
  m <- length(X)

  if (is.null(cc)) {

      if ((method == "MLE")||(method == "MLE+MOM")) {
        method1 <- "CSS-ML"
      } else {
        method1 <- "CSS"
      }
    
    if (case == 'U') {
      model <- arima(X, order = neworder, method = method1)
      
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
  
      cc <- getCC.ARMA(fap0 = fap0, interval = interval, m, order = order, phi.vec = phi.vec, theta.vec = theta.vec, case = case,
        method = method, nsim.coefs = nsim.coefs, nsim.process = nsim.process, burn.in = burn.in, sim.type = sim.type, verbose = verbose)

    } else if (case == 'K') {
      
      cc <- getCC.ARMA(fap0 = fap0, interval = interval, m, order = order, phi.vec = phi.vec, theta.vec = theta.vec, case = case,
        method = method, nsim.coefs = nsim.coefs, nsim.process = nsim.process, burn.in = burn.in, sim.type = sim.type, verbose = verbose)

    }
    
    
  }
  

  if (case == "U") {
     mu <- mean(X)
     gamma <- sd(X)
      
  } else if (case == "K") {
    mu <- mu0
    gamma <- sigma0
    
  }
  
  if (standardize){
    CS <- (X - mu) / gamma
  
    LCL <- -cc
    UCL <- cc
  }else{
    CS <- X
  
    LCL <- -cc * gamma + mu
    UCL <- cc * gamma + mu
    
    if (transform == "boxcox") {
      CS <- forecast::InvBoxCox(CS, lambda1)
      mu <- forecast::InvBoxCox(mu, lambda1)
      LCL <- forecast::InvBoxCox(LCL, lambda1)
      UCL <- forecast::InvBoxCox(UCL, lambda1)
    } else if (transform == "yeojohnson") {
      CS <- VGAM::yeo.johnson(CS, lambda1, inverse = TRUE)
      mu <- VGAM::yeo.johnson(mu, lambda1, inverse = TRUE)
      LCL <- VGAM::yeo.johnson(LCL, lambda1, inverse = TRUE)
      UCL <- VGAM::yeo.johnson(UCL, lambda1, inverse = TRUE)
    }
  }
  
  
  if (anyNA(CS) | is.na(mu) | is.na(LCL)| is.na(UCL)) {
    if (transform != 'none') {
      warning("The back transformation is erroneous using the given lambda.  Please try other lambda.")
    }
  }
  
  
  if (plot.option == TRUE) {

    main.text <- paste('Phase I Individual Chart')

    plot(c(1, m), c(min(LCL, CS, na.rm = TRUE), max(UCL, CS, na.rm = TRUE)), xaxt = "n", xlab = 'Observation', ylab = 'Charting Statistic', type = 'n', main = main.text)

    axis(side = 1, at = 1:m)

    points(1:m, CS, type = 'o')
    points(c(-1, m + 2), c(LCL, LCL), type = 'l', col = 'red')
    points(c(-1, m + 2), c(UCL, UCL), type = 'l', col = 'red')
    if (standardize) {
      points(c(-1, m + 2), c(0, 0), type = 'l', col = 'blue')
    } else {
      points(c(-1, m + 2), c(mu, mu), type = 'l', col = 'blue')
    }
    
    text(round(m * 0.8), UCL, paste('UCL = ', round(UCL, 4)), pos = 1)
    text(round(m * 0.8), LCL, paste('LCL = ', round(LCL, 4)), pos = 3)

  }

  res <- list(CL = mu, gamma = gamma, cc = cc, order = order, phi.vec = phi.vec, theta.vec = theta.vec, 
              LCL = LCL, UCL = UCL, CS = CS, transform = transform, lambda = lambda1, standardize = standardize)

  return(invisible(res))
}
