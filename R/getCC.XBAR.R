getcc.xbar <- function(fap0,
                       m,
                       var.est = c("S", "MR"),
                       ub.cons = 1,
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
      fap0 = fap0, interval = interval, m = m, est = var.est,
      ub.cons = ub.cons, nsim = nsim, verbose = verbose
    )
  } else if (method == "BA") {
    getcc.ba(fap0 = fap0, m = m, nu = nu, ub.cons = ub.cons, lambda = lambda)
  } else {
    stop('Unexpected variance estimator. The method must be in c("exact", "BA"). The program will stop.')
  }
}

getcc <- function(fap0,
                  m,
                  var.est = c("S", "MR"),
                  ub.cons = 1,
                  method = c("exact", "BA"),
                  interval = c(1, 5),
                  nsim = 10000,
                  nu = m - 1,
                  lambda = 1,
                  verbose = FALSE) {
  warning("This function has been renamed.  Please use getcc.xbar instead.")

  getcc.xbar(
    fap0 = fap0,
    m = m,
    var.est = var.est,
    ub.cons = ub.cons,
    method = method,
    interval = interval,
    nsim = nsim,
    nu = nu,
    lambda = lambda,
    verbose = verbose
  )
}
