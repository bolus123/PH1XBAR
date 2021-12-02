sigma.mat.f <- function(n) {
  sigmaMat <- diag(n)
  sigmaMat[which(sigmaMat == 0)] <- -1 / (n - 1)
  return(sigmaMat)
}

integrand.s <- function(Z, n, c, ub.cons = 1) {
  Z <- as.vector(Z)

  Numerator <- Z^2
  Denominator <- sum(Numerator) / (n - 1)

  NSE <- Numerator / Denominator < (c / ub.cons)^2

  idicator <- ifelse(sum(NSE) == n, 1, 0)

  return(idicator)
}

integrand.mr.bar <- function(Z, n, c, ub.cons = 1) {
  Z <- as.vector(Z)

  Numerator <- Z
  Denominator <- mean(abs(diff(Z)))

  NSE <- (Numerator / Denominator)^2 < (c / ub.cons)^2

  idicator <- ifelse(sum(NSE) == n, 1, 0)

  return(idicator)
}

integral.mc <- function(n, c, est = c("S", "MR"), ub.cons = 1, nsim = 10000) {
  sigmaMat <- sigma.mat.f(n)

  res <- lapply(X = 1:nsim, function(X) {
    Z <- rmvnorm(1, sigma = sigmaMat)
    if (est == "S") {
      integrand.s(Z, n, c, ub.cons)
    } else if (est == "MR") {
      integrand.mr.bar(Z, n, c, ub.cons)
    }
  })

  1 - mean(unlist(res))
}

getCC.exact.n <- function(fap0, interval = c(1, 5), n, est = c("S", "MR"), ub.cons = 1,
                          nsim = 10000, verbose = FALSE) {
  root.finding <- function(c, fap0, n, est = c("S", "MR"), ub.cons = 1,
                           nsim = 10000) {
    fapin <- integral.mc(n, c, est, ub.cons, nsim)
    if (verbose) {
      cat("fapin:", fapin, " and cc:", c, "\n")
    }
    fap0 - fapin
  }

  est <- est[1]

  uniroot(
    f = root.finding, interval = interval, fap0 = fap0, n = n, est = est,
    ub.cons = ub.cons, nsim = nsim
  )$root
}

getCC.exact <- function(fap0, interval = c(1, 5), m, est = c("S", "MR"), ub.cons = 1,
                        nsim = 10000, verbose = FALSE) {
  est <- est[1]

  getCC.exact.n(fap0, interval, m,
    est = est, ub.cons = ub.cons,
    nsim = nsim, verbose = verbose
  )
}
