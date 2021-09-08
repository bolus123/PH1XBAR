sigmaMatF <- function(n) {
  
  sigmaMat <- diag(n)
  sigmaMat[which(sigmaMat == 0)] <- - 1 / (n - 1)
  return(sigmaMat)
  
}

integrandS <- function(Z, n, c, ubCons = 1) {
  
  Z <- as.vector(Z)
  
  Numerator <- Z ^ 2
  Denominator <- sum(Numerator) / (n - 1)
  
  NSE <- Numerator / Denominator < (c / ubCons) ^ 2
  
  idicator <- ifelse(sum(NSE) == n, 1, 0)

  return(idicator)
  
}

integrandMRBar <- function(Z, n, c, ubCons = 1) {
  
  Z <- as.vector(Z)
  
  Numerator <- Z
  Denominator <- mean(abs(diff(Z)))
  
  NSE <- (Numerator / Denominator) ^ 2 < (c / ubCons) ^ 2
  
  idicator <- ifelse(sum(NSE) == n, 1, 0)
  
  return(idicator)
  
}

integralMC <- function(n, c, est = c('S', 'MR'), ubCons = 1, nsim = 10000) {
  
  sigmaMat <- sigmaMatF(n)
  
  res <- lapply(X = 1:nsim, function(X){
    Z <- rmvnorm(1, sigma = sigmaMat)
    if (est == 'S') {
      integrandS(Z, n, c, ubCons)
    } else if (est == 'MR') {
      integrandMRBar(Z, n, c, ubCons)
    }
  })
  
  1 - mean(unlist(res))
}

getCC.exact.n <- function(FAP0, interval = c(1, 5), n, est = c('S', 'MR'), ubCons = 1, 
                         nsim = 10000) {
  
  root.finding <- function(c, FAP0, n, est = c('S', 'MR'), ubCons = 1, 
                           nsim = 10000){
    FAPin <- integralMC(n, c, est, ubCons, nsim)
    cat('FAPin:', FAPin, ' and cc:', c, '\n')
    FAP0 - FAPin
  }
  
  est <- est[1]
  
  uniroot(f = root.finding, interval = interval, FAP0 = FAP0, n = n, est = est,
          ubCons = ubCons, nsim = nsim)$root
  
}

getCC.exact <- function(FAP0, interval = c(1, 5), m, est = c('S', 'MR'), ubCons = 1, 
                         nsim = 10000) {
		
	est <- est[1]
		
	getCC.exact.n(FAP0, interval, m, est = est, ubCons = ubCons, 
                         nsim = nsim)
		
}
