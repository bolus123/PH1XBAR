getCC <- function(
      FAP0
			,m
			,var.est = c('S', 'MR')
			,ubCons = 1 
			,method = c('exact', 'BA')
			,interval = c(1, 5)
			,nsim = 10000
			,nu = m - 1 
			,lambda = 1

  ){
   
  var.est <- var.est[1]
  method <- method[1]
  
  if (method == 'exact') {
      
    getCC.exact(FAP0 = FAP0, interval = interval, m = m, est = var.est, 
                ubCons = ubCons, nsim = nsim)
  
  } else if (method == 'BA') {
  
    getCC.BA(FAP0 = FAP0, m = m, nu = nu, ubCons = ubCons, lambda = lambda)
  
  } else {
    
    stop('Unexpected variance estimator. The method must be in c("exact", "BA"). The program will stop.')
  
  }

}



getCC.ARMA <- function(FAP0 = 0.1, interval = c(1, 4), n = 50, mvVec = NULL, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, case = 'U',
              method = 'Method 3', nsimCoefs = 100, nsimProcess = 1000, burnIn = 50, simType = 'Matrix', logliktol = 1e-2, seed = 12345) {


  if (case == 'K') {
    if (simType == 'Matrix') {
      sigMat <- SigmaMat(n, order = order, phiVec = phiVec, thetaVec = thetaVec, sigma2 = 1, burnIn = burnIn)
      out <- qmvnorm(1 - FAP0, tail = 'both.tails', sigma = sigMat$SigmaMat / sigMat$gamma0)$quantile
    } else {
      if (is.null(mvVec)) {
        out <- getCCPH1ARMASim(FAP0, interval, n, order, phiVec = phiVec, thetaVec = thetaVec, case = case, method = method,
                               nsim = nsimProcess, burnIn = burnIn, simType = simType, seed = seed)
      } else {
        out <- getCCPH1ARMASimMissingValue(FAP0, interval, n, mvVec, order, phiVec = phiVec, thetaVec = thetaVec, case = case, method = method,
                               nsim = nsimProcess, burnIn = burnIn, simType = simType, logliktol = logliktol)
      }

    }
  } else if (case == 'U') {

    if (is.null(mvVec)) {
      CoefDist <- simCoefDist(n, order, phiVec, thetaVec, method, nsim = nsimCoefs, burnIn = burnIn, simType = simType)

      out <- getCCPH1ARMASim(FAP0, interval, n, order, phiVec = CoefDist$phiVec, thetaVec = CoefDist$thetaVec, case = case, method = method,
                             nsim = nsimProcess, burnIn = burnIn, simType = simType, seed = seed)
    } else {

      CoefDist <- simCoefDistMissingValue(n, mvVec, order, phiVec, thetaVec, method, nsim = nsimCoefs, burnIn = burnIn, tol = logliktol, simType = simType)

      out <- getCCPH1ARMASimMissingValue(FAP0, interval, n, mvVec, order, phiVec = CoefDist$phiVec, thetaVec = CoefDist$thetaVec, case = case, method = method,
                             nsim = nsimProcess, burnIn = burnIn, simType = simType, logliktol = logliktol)
    }

  }

  out

}