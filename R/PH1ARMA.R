
PH1ARMA <- function(X, cc = NULL, FAP0 = 0.1, order = NULL, plot.option = TRUE, interval = c(1, 4),
                    case = 'U', method = 'Method 3', nsimCoefs = 100, nsimProcess = 1000, burnIn = 50, 
                    simType = 'Matrix', logliktol = 1e-2, verbose = FALSE) {

  Y <- X

  n <- length(Y)

  if (is.null(cc)) {

    if (is.null(order)) {

      if (method == 'Method 1' | method == 'Method 3') {
        model <- auto.arima(Y, method = 'CSS-ML')
      } else if (method == 'Method 2') {
        model <- auto.arima(Y, method = 'CSS')
      }

      order <- rep(0, 3)
      order[1] <- length(model$model$phi)
      order[2] <- length(model$model$Delta)
      order[3] <- length(model$model$theta)

    } else {

      model <- arima(Y, order = order, method = method)

    }

    if (length(model$model$phi) > 0) {
      phiVec <- model$model$phi
    } else {
      phiVec <- NULL
    }

    if (length(model$model$theta) > 0) {
      thetaVec <- model$model$theta
    } else {
      thetaVec <- NULL
    }

    cc <- getCC.ARMA(FAP0 = FAP0, interval = interval, n, order = order, phiVec = phiVec, thetaVec = thetaVec, case = case,
      method = method, nsimCoefs = nsimCoefs, nsimProcess = nsimProcess, burnIn = burnIn, simType = simType, verbose = verbose)

  }

  if (order[2] > 0) {

    Y <- diff(Y, differences = order[2])

  }

  mu <- mean(Y)
  gamma <- sd(Y)

  stdX <- (Y - mu) / gamma

  LCL <- -cc
  UCL <- cc

  if (plot.option == TRUE) {

    main.text <- paste('Phase I Individual Chart for FAP0 =', FAP0, 'with an ARMA model')

    plot(c(1, n), c(min(LCL, stdX), max(UCL, stdX)), xaxt = "n", xlab = 'Observation', ylab = 'Charting Statistic', type = 'n', main = main.text)

    axis(side = 1, at = 1:n)

    points(1:n, stdX, type = 'o')
    points(c(-1, n + 2), c(LCL, LCL), type = 'l', col = 'red')
    points(c(-1, n + 2), c(UCL, UCL), type = 'l', col = 'red')
    points(c(-1, n + 2), c(mu, mu), type = 'l', col = 'blue')
    text(round(n * 0.8), UCL, paste('UCL = ', round(UCL, 4)), pos = 1)
    text(round(n * 0.8), LCL, paste('LCL = ', round(LCL, 4)), pos = 3)

  }

  list(CL = mu, gamma0 = gamma, cc = cc, order = order, phiVec = phiVec, thetaVec = thetaVec, LCL = LCL, UCL = UCL, CS = stdX)


}
