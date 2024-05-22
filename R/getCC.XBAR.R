#' Random Flexible Level Shift Model
#' 
#' get Phase I corrected charting constant
#' 
#' @aliases getCC.XBAR
#' @aliases getCC
#' 
#' @param m nominal false Alarm Probabilty in Phase 1
#' @param fap0 number of subgroups
#' @param var.est 'S' - use mean-square-based estimator, 'MR' - use moving-range-based estimator
#' @param ub.cons unbiasing constant 
#' @param method 'exact' - calculate results using the exact method, 'BA' - calculate results using the Bonfferoni approximation 
#' @param interval searching range of charting constants for the exact method 
#' @param nsim number of simulation for the exact method 
#' @param nu degrees of freedom for the Bonfferoni approximation
#' @param lambda constant for the Bonfferoni approximation 
#' @param verbose print diagnostic information about fap0 and the charting constant during the simulations for the exact method 
#' @return Object type double. The corrected charting constant. 
#' 
#' @export
#' @examples
#' set.seed(12345)
#' 
#' # Calculate the charting constant using 10 simulations and mean-square-based estimator
#' getCC.XBAR(fap0=0.05, m=20, nsim=10, var.est='S', verbose = TRUE)
#' 
#' # Calculate the charting constant using 10 simulations and moving-range-based estimator
#' getCC.XBAR(fap0=0.05, m=20, nsim=10, var.est='MR', verbose = TRUE)
#' 
#' 
getCC.XBAR <- function(m,
                       fap0 = 0.05,
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
    getCC.exact(
      fap0 = fap0, interval = interval, m = m, est = var.est,
      ub.cons = ub.cons, nsim = nsim, verbose = verbose
    )
  } else if (method == "BA") {
    getCC.ba(fap0 = fap0, m = m, nu = nu, ub.cons = ub.cons, lambda = lambda)
  } else {
    stop('Unexpected variance estimator. The method must be in c("exact", "BA"). The program will stop.')
  }
}

#' @export
getCC <- function(m,
                  fap0 = 0.05,
                  var.est = c("S", "MR"),
                  ub.cons = 1,
                  method = c("exact", "BA"),
                  interval = c(1, 5),
                  nsim = 10000,
                  nu = m - 1,
                  lambda = 1,
                  verbose = FALSE) {
  warning("This function has been renamed.  Please use getCC.XBAR instead.")

  getCC.XBAR(
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
