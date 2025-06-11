##############################################
# Matrix form
##############################################
invert.q <- function(coef) {
  out <- 1

  if (all(abs(coef) < 1)) {
    minmod <- min(Mod(polyroot(c(1, coef))))

    if (minmod <= 1) {
      out <- 0
    }
  } else {
    out <- 0
  }

  return(out)
}

pars.mat <- function(n, parsVec, norder = 1) {
  Check <- invert.q(parsVec)
  if (Check == 0) {
    NULL
  } else {
    Mat <- diag(n)
    for (i in 1:norder) {
      Mat <- Mat + Diag(rep(parsVec[i], n - i), k = -i)
    }
    Mat
  }
}



sigma.mat <- function(n, order = c(1, 0, 0), phi.vec = 0.5, theta.vec = NULL, sigma2 = 1, burn.in = 50) {
  if (order[1] == 0) {
    phiMat <- diag(n + burn.in)
  } else {
    phiMat <- pars.mat(n + burn.in, -phi.vec, norder = order[1])
  }

  if (order[3] == 0) {
    thetaMat <- diag(n + burn.in)
  } else {
    thetaMat <- pars.mat(n + burn.in, theta.vec, norder = order[3])
  }

  out <- solve(phiMat) %*% thetaMat %*% t(thetaMat) %*% t(solve(phiMat)) * sigma2
  

  gamma0 <- out[dim(out)[1], dim(out)[2]]

  if (burn.in > 0) {
    out <- out[-c(1:burn.in), -c(1:burn.in)]
  }

  list(sigma.mat = out, sqrtsigma.mat = chol(out), gamma0 = gamma0)
}

##############################################
# process simulation
##############################################

# Simulate innovations
sim.innov <- function(n, sigma2 = 1, XSim = "norm", XPars = c(0, 1)) {
  if (XSim == "norm") {
    me <- XPars[1]
    std <- sqrt(XPars[2])

    pars <- c(me, std)

    rgen <- rnorm
  } else if (XSim == "t") {
    me <- 0
    std <- sqrt(XPars[1] / (XPars[1] - 2))

    pars <- c(XPars[1])

    rgen <- rt
  } else if (XSim == "gamma") {
    me <- XPars[1] * XPars[2]
    std <- sqrt(XPars[1] * XPars[2]^2)

    pars <- c(XPars[1], 1 / XPars[2])

    rgen <- rgamma
  } else if (XSim == "beta") {
    me <- XPars[1] / (XPars[1] + XPars[2])
    std <- sqrt(XPars[1] * XPars[2] / ((XPars[1] + XPars[2])^2) / (XPars[1] + XPars[2] + 1))

    pars <- c(XPars[1], XPars[2])

    rgen <- rbeta
  }

  if (XSim != "t") {
    out <- rgen(n, pars[1], pars[2])
  } else {
    out <- rgen(n, pars[1])
  }

  (out - me) / std * sqrt(sigma2)
}



sim.ARMA.process <- function(n, order = c(1, 0, 0), phi.vec = 0.5, theta.vec = NULL, sigma2 = 1, innovDist = "norm", innovPars = c(0, 1), burn.in = 50,
                           sim.type = "Matrix", SigMat = sigma.mat(n = n, order = order, phi.vec = phi.vec, theta.vec = theta.vec, sigma2 = sigma2, burn.in = burn.in)) {
  if (sim.type == "Matrix") {
    as.vector(SigMat$sqrtsigma.mat %*% matrix(sim.innov(n, sigma2 = 1, XSim = innovDist, XPars = innovPars), ncol = 1))
  } else if (sim.type == "Recursive") {
    innov <- sim.innov(n + burn.in + order[1] + order[3], sigma2 = sigma2, XSim = innovDist, XPars = innovPars)

    n.start <- order[1] + order[3]

    if (burn.in > 0) {
      arima.sim(list(order = order, ar = phi.vec, ma = theta.vec),
        n = n + burn.in, innov = innov[-c(1:(n.start))], n.start = n.start, start.innov = innov[c(1:(n.start))]
      )[(burn.in + 1):(n + burn.in)]
    } else {
      arima.sim(list(order = order, ar = phi.vec, ma = theta.vec),
        n = n + burn.in, innov = innov[-c(1:(n.start))], n.start = n.start, start.innov = innov[c(1:(n.start))]
      )
    }
  }
}

sim.coef.dist <- function(n, order = c(1, 0, 0), phi.vec = 0.5, theta.vec = NULL, method = "CSS",
                        nsim = 100, burn.in = 50, sim.type = "Matrix",
                        SigMat = sigma.mat(n, order, phi.vec, theta.vec, sigma2 = 1, burn.in = burn.in)) {
  outAR <- matrix(NA, nrow = nsim, ncol = order[1])
  outMA <- matrix(NA, nrow = nsim, ncol = order[3])
  outMean <- rep(NA, nsim)
  outGamma <- rep(NA, nsim)

  for (i in 1:nsim) {
    flg <- 1

    while (flg == 1) {
      sim <- sim.ARMA.process(n, order, phi.vec, theta.vec,
        sigma2 = 1, innovDist = "norm", innovPars = c(0, 1), burn.in = burn.in,
        sim.type = sim.type, SigMat = SigMat
      )

      ##if (method == "MLE" | method == "MLE+MOM") {
      ##  model <- try(arima(sim, order = order, method = "CSS-ML"), silent = TRUE)
      ##} else if (method == "CSS") {
      ##  model <- try(arima(sim, order = order, method = "CSS"), silent = TRUE)
      ##}
      
      model <- try(arima(sim, order = order, method = method), silent = TRUE)

      check1 <- 1
      check2 <- 1

      if (class(model)[1] != "try-error") {
        if (order[1] > 0) {
          outAR[i, ] <- model$coef[1:order[1]]
          check1 <- invert.q(outAR[i, ]) == 1
          # cat('AR:', outAR[i, ], '\n')
          # cat('check1:', check1, '\n')
        } else {
          outAR[i, ] <- rep(0, order[1])
        }

        if (order[3] > 0) {
          outMA[i, ] <- model$coef[(order[1] + 1):(order[1] + order[3])]
          check2 <- invert.q(outMA[i, ]) == 1
          # cat('MA:', outMA[i, ], '\n')
          # cat('check2:', check2, '\n')
        } else {
          outMA[i, ] <- rep(0, order[3])
        }


        if (check1 & check2) {
          outMean[i] <- mean(sim)
          outGamma[i] <- var(sim)

          flg <- 0
        }
      }
    }
  }

  list(phi.vec = outAR, theta.vec = outMA)
}




fap.PH1ARMA <- function(cc = 3, n = 50, order = c(1, 0, 0), phi.vec = 0.5, theta.vec = NULL, case = "U", method = "MLE+MOM",
                       nsim = 1000, burn.in = 50, sim.type = "Matrix") {
  if (sim.type == "Matrix") {
    sigMat1 <- sigma.mat(n, order, phi.vec, theta.vec, sigma2 = 1, burn.in = burn.in)
    gamma0 <- sigMat1$gamma0
  } else if (sim.type == "Recursive") {
    if (case == "K") {
      sigMat1 <- sigma.mat(100, order, phi.vec, theta.vec, sigma2 = 1, burn.in = 50)
      gamma0 <- sigMat1$gamma0
    } else {
      sigMat1 <- NULL
    }
  }

  out <- lapply(1:nsim, function(X) {
    flg <- 1
    while (flg == 1) {
      sim <- sim.ARMA.process(n, order, phi.vec, theta.vec,
        sigma2 = 1,
        innovDist = "norm", innovPars = c(0, 1), burn.in = burn.in, sim.type = sim.type, SigMat = sigMat1
      )

      if (case == "U") {
        if (method == "CSS" | method == "MLE+MOM") {
          sim <- (sim - mean(sim)) / sd(sim)
          flg <- 0
          out1 <- sum(-cc <= sim & sim <= cc) != n
          return(out1)
        } else if (method == "MLE") {
          m1 <- try(arima(sim, order = order, method = "CSS-ML"), silent = TRUE)

          if (class(m1)[1] != "try-error") {
            if (!is.null(phi.vec)) {
              phi.vecNew <- m1$coef[1:order[1]]
            } else {
              phi.vecNew <- NULL
            }

            if (!is.null(theta.vec)) {
              theta.vecNew <- m1$coef[(order[1] + 1):(order[1] + order[3])]
            } else {
              theta.vecNew <- NULL
            }

            Intercept <- m1$coef[order[1] + order[3] + 1]
            mu0 <- Intercept * (1 - sum(phi.vecNew))
            sigma2 <- m1$sigma2
            gamma0 <- sigma.mat(100, order = order, phi.vec = phi.vec, theta.vec = theta.vec, sigma2 = sigma2, burn.in = 50)$gamma0

            sim <- (sim - mu0) / sqrt(gamma0)
            flg0 <- 0
            out1 <- sum(-cc <= sim & sim <= cc) != n
            return(out1)
          }
        }
      } else if (case == "K") {
        sim <- (sim) / sqrt(gamma0)
        flg <- 0
        out1 <- sum(-cc <= sim & sim <= cc) != n
        return(out1)
      }
    }
  })

  mean(unlist(out))
}


getCC.PH1ARMA.sim <- function(fap0 = 0.1, interval = c(1, 4), n = 50, order = c(1, 0, 0), phi.vec = 0.5, theta.vec = NULL, case = "U", method = "MLE+MOM",
                            nsim = 1000, burn.in = 50, sim.type = "Matrix", verbose = FALSE) {
  root.finding <- function(fap0, cc, n, order, phi.vec, theta.vec, case, method, nsim1, nsim2, burn.in, sim.type) {
    if (nsim1 > 0) {
      fapin <- lapply(1:nsim1, function(X) {
        phi.vecTmp <- phi.vec[X, ]
        if (length(phi.vecTmp) == 0) phi.vecTmp <- rep(0, order[1])

        theta.vecTmp <- theta.vec[X, ]
        if (length(theta.vecTmp) == 0) theta.vecTmp <- rep(0, order[3])

        # cat('AR:', phi.vecTmp, '\n')
        # cat('checkAR:', all(abs(outAR[i, ]) < 1), '\n')

        # cat('MA:', theta.vecTmp, '\n')
        # cat('checkMA:', all(abs(outMA[i, ]) < 1), '\n')

        fap.PH1ARMA(
          cc = cc, n = n, order = order, phi.vec = phi.vecTmp, theta.vec = theta.vecTmp,
          case = case, method = method, nsim = nsim2, burn.in = burn.in, sim.type = sim.type
        )
      })
      fapin <- mean(unlist(fapin))
    } else {
      fapin <- fap.PH1ARMA(
        cc = cc, n = n, order = order, phi.vec = phi.vec, theta.vec = theta.vec,
        case = case, method = method, nsim = nsim2, burn.in = burn.in, sim.type = sim.type
      )
    }

    if (verbose) {
      cat("fapin:", fapin, " and cc:", cc, "\n")
    }
    fap0 - fapin
  }

  if (is.matrix(phi.vec) | is.matrix(theta.vec)) {
    nsim1 <- max(dim(phi.vec)[1], dim(theta.vec)[1])
  } else {
    nsim1 <- 0
  }

  ## cat('phi.vec:', phi.vec, 'theta.vec:', theta.vec, sep = ', ')

  uniroot(root.finding, interval,
    fap0 = fap0, n = n, order = order, phi.vec = phi.vec, theta.vec = theta.vec, case = case, method = method,
    nsim1 = nsim1, nsim2 = nsim, burn.in = burn.in, sim.type = sim.type
  )$root
}


#' get Phase I corrected charting constant with an ARMA model
#' 
#' @param fap0 nominal false Alarm Probabilty in Phase 1 
#' @param interval searching range of charting constants for the exact method 
#' @param n number of observations 
#' @param order order for ARMA model 
#' @param phi.vec given vectors of autoregressive parameters for ARMA models 
#' @param theta.vec given vectors of moving-average parameters for ARMA models 
#' @param case known or unknown case.  When case = 'U', the parameters are unknown and the charting constant is calculated based on a bootstrapping method.  When case = 'K', the parameters are known and the charting constant is found using the quantile function of multivariate normal distribution
#' @param phi vector of autoregressive coefficient(s).  When case = 'K', it must be provided. The length must be the same as the first value in the order.  It needs to be NULL if no autoregressor presents
#' @param theta vector of moving-average coefficient(s).  When case = 'K', it must be provided. The length must be the same as the third value in the order.  It needs to be NULL if no moving average presents
#' @param method estimation method for the control chart. When method = 'Method 3' is maximum likehood estimations plus method of moments. Other options are 'Method 1' which is pure MLE and 'Method 2' which is pure CSS. 
#' @param nsim.coefs number of simulation for coeficients.  It is functional when double.sim = TRUE. 
#' @param nsim.process number of simulation for ARMA processes 
#' @param burn.in number of burn-ins.  When burn.in = 0, the ECM gets involved.  When burn.in is large enough, the ACM gets involved. 
#' @param sim.type type of simulation.  When sim.type = 'Matrix', the simulation is generated using matrix computation.  When sim.type = 'Recursive', the simulation is based on a recursion. 
#' @param verbose print diagnostic information about fap0 and the charting constant during the simulations for the exact method 
#' @return Object type double. The corrected charting constant.
#' 
#' 
#' @export
#' 
#' @examples
#' # load the data in the package as an example
#' set.seed(12345)
#' 
#' # Calculate the charting constant using fap0 of 0.05, and 50 observations
#' getCC.ARMA(fap0=0.05, n=50, nsim.coefs=10, nsim.process=10)
#' 
getCC.ARMA <- function(fap0 = 0.05,
                       interval = c(1, 4),
                       n = 50,
                       order = c(1, 0, 0),
                       phi.vec = 0.5,
                       theta.vec = NULL,
                       case = "U",
                       method = "MLE+MOM",
                       nsim.coefs = 100,
                       nsim.process = 1000,
                       burn.in = 50,
                       sim.type = "Matrix",
                       verbose = FALSE) {
  if (case == "K") {
    if (sim.type == "Matrix") {
      sigMat <- sigma.mat(n, order = order, phi.vec = phi.vec, theta.vec = theta.vec, sigma2 = 1, burn.in = burn.in)
      out <- qmvnorm(1 - fap0, tail = "both.tails", sigma = sigMat$sigma.mat / sigMat$gamma0)$quantile
    } else {
      out <- getCC.PH1ARMA.sim(fap0, interval, n, order,
        phi.vec = phi.vec, theta.vec = theta.vec, case = case, method = method,
        nsim = nsim.process, burn.in = burn.in, sim.type = sim.type, verbose = verbose
      )
    }
  } else if (case == "U") {
    CoefDist <- sim.coef.dist(n, order, phi.vec, theta.vec, method, nsim = nsim.coefs, burn.in = burn.in, sim.type = sim.type)

    out <- getCC.PH1ARMA.sim(fap0, interval, n, order,
      phi.vec = CoefDist$phi.vec, theta.vec = CoefDist$theta.vec, case = case, method = method,
      nsim = nsim.process, burn.in = burn.in, sim.type = sim.type, verbose = verbose
    )
  }

  out
}
