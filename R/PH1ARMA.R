

ph1arma <- function(X, cc = NULL, fap0 = 0.1, order = NULL, plot.option = TRUE, interval = c(1, 4),
                    case = 'U', method = 'Method 3', nsim.coefs = 100, nsim.process = 1000, burn.in = 50, 
                    sim.type = 'Matrix', logliktol = 1e-2, verbose = FALSE) {
  
  if (!is.vector(X)) {
	if (dim(X)[1] == 1 | dim(X)[2] == 1) {
		X <- as.vector(X)
	} else {
		stop('X is not a vector, or a m x 1 or 1 x m matrix.')
	}
  } 

  n <- length(X)

  if (is.null(cc)) {

    if (is.null(order)) {

      if (method == 'Method 1' | method == 'Method 3') {
        model <- auto.arima(X, method="CSS-ML", max.p=1, max.q=0, max.d=1)
      } else if (method == 'Method 2') {
        model <- auto.arima(X, method="CSS", max.p=1, max.q=0, max.d=1)
      }

      order <- rep(0, 3)
      order[1] <- length(model$model$phi)
      order[2] <- length(model$model$Delta)
      order[3] <- length(model$model$theta)

	  if (sum(model$model$phi) == 0) {
		order[1] <- 0
	  }
	  
	  if (sum(model$model$theta) == 0) {
		order[3] <- 0
	  }

    } else {

      model <- arima(X, order = order, method = method)

    }

    if (length(model$model$phi) > 0) {
      if (any(model$model$phi > 1.0)){
        model$model$phi[model$model$phi >= 1.0] = 0.99
        model$model$phi[model$model$phi <= -1.0] = -0.99
      }
      phi.vec <- model$model$phi
    } else {
      phi.vec <- NULL
    }

    if (length(model$model$theta) > 0) {
      if (any(model$model$theta > 1.0)){
        model$model$theta[model$model$theta >= 1.0] = 0.99
        model$model$theta[model$model$theta <= -1.0] = -0.99
      }
      theta.vec <- model$model$theta
    } else {
      theta.vec <- NULL
    }

    cc <- getcc.arma(fap0 = fap0, interval = interval, n, order = order, phi.vec = phi.vec, theta.vec = theta.vec, case = case,
      method = method, nsim.coefs = nsim.coefs, nsim.process = nsim.process, burn.in = burn.in, sim.type = sim.type, verbose = verbose)

  }

  if (order[2] > 0) {

    X <- diff(X, differences = order[2])

  }

  mu <- mean(X)
  gamma <- sd(X)

  stdX <- (X - mu) / gamma

  LCL <- -cc
  UCL <- cc

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
