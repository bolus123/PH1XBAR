

PH1ARMA <- function(X, cc = NULL, fap0 = 0.05, order = c(1, 0, 0), plot.option = TRUE, interval = c(1, 4),
                    case = 'U', phiVec = NULL, thetaVec = NULL, mu0 = NULL, sigma0 = NULL, method = 'MLE+MOM', nsim.coefs = 100, nsim.process = 1000, burn.in = 50, 
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
      phi.vec <- phiVec
      theta.vec <- thetaVec
      
      cc <- getCC.ARMA(fap0 = fap0, interval = interval, n, order = order, phi.vec = phi, theta.vec = theta, case = case,
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
