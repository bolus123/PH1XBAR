sigma.mat.f <- function(n) {
  sigmaMat <- diag(n)
  sigmaMat[which(sigmaMat == 0)] <- -1 / (n - 1)
  return(sigmaMat)
}

integrand.s <- function(Z, n, c, ubCons = 1) {
  Z <- as.vector(Z)

  Numerator <- Z^2
  Denominator <- sum(Numerator) / (n - 1)

  NSE <- Numerator / Denominator < (c / ubCons)^2

  idicator <- ifelse(sum(NSE) == n, 1, 0)

  return(idicator)
}

integrand.mr.bar <- function(Z, n, c, ubCons = 1) {
  Z <- as.vector(Z)

  Numerator <- Z
  Denominator <- mean(abs(diff(Z)))

  NSE <- (Numerator / Denominator)^2 < (c / ubCons)^2

  idicator <- ifelse(sum(NSE) == n, 1, 0)

  return(idicator)
}

integral.mc <- function(n, c, est = c("S", "MR"), ubCons = 1, nsim = 10000) {
  sigmaMat <- sigma.mat.f(n)

  res <- lapply(X = 1:nsim, function(X) {
    Z <- rmvnorm(1, sigma = sigmaMat)
    if (est == "S") {
      integrand.s(Z, n, c, ubCons)
    } else if (est == "MR") {
      integrand.mr.bar(Z, n, c, ubCons)
    }
  })

  1 - mean(unlist(res))
}

getcc.exact.n <- function(FAP0, interval = c(1, 5), n, est = c("S", "MR"), ubCons = 1,
                          nsim = 10000, verbose = FALSE) {
  root.finding <- function(c, FAP0, n, est = c("S", "MR"), ubCons = 1,
                           nsim = 10000) {
    FAPin <- integral.mc(n, c, est, ubCons, nsim)
    if (verbose) {
      cat("FAPin:", FAPin, " and cc:", c, "\n")
    }
    FAP0 - FAPin
  }

  est <- est[1]

  uniroot(
    f = root.finding, interval = interval, FAP0 = FAP0, n = n, est = est,
    ubCons = ubCons, nsim = nsim
  )$root
}

getcc.exact <- function(FAP0, interval = c(1, 5), m, est = c("S", "MR"), ubCons = 1,
                        nsim = 10000, verbose = FALSE) {
  est <- est[1]

  getcc.exact.n(FAP0, interval, m,
    est = est, ubCons = ubCons,
    nsim = nsim, verbose = verbose
  )
}
