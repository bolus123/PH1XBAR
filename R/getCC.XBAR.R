getcc.xbar <- function(FAP0,
                       m,
                       var.est = c("S", "MR"),
                       ubCons = 1,
                       method = c("exact", "BA"),
                       interval = c(1, 5),
                       nsim = 10000,
                       nu = m - 1,
                       lambda = 1,
                       verbose = FALSE) {
  var.est <- var.est[1]
  method <- method[1]

  if (method == "exact") {
    getcc.exact(
      FAP0 = FAP0, interval = interval, m = m, est = var.est,
      ubCons = ubCons, nsim = nsim, verbose = verbose
    )
  } else if (method == "BA") {
    getcc.ba(FAP0 = FAP0, m = m, nu = nu, ubCons = ubCons, lambda = lambda)
  } else {
    stop('Unexpected variance estimator. The method must be in c("exact", "BA"). The program will stop.')
  }
}

getcc <- function(FAP0,
                  m,
                  var.est = c("S", "MR"),
                  ubCons = 1,
                  method = c("exact", "BA"),
                  interval = c(1, 5),
                  nsim = 10000,
                  nu = m - 1,
                  lambda = 1,
                  verbose = FALSE) {
  warning("This function has been renamed.  Please use getcc.xbar instead.")

  getcc.xbar(
    FAP0 = FAP0,
    m = m,
    var.est = var.est,
    ubCons = ubCons,
    method = method,
    interval = interval,
    nsim = nsim,
    nu = nu,
    lambda = lambda,
    verbose = verbose
  )
}
