#' Random Flexible Level Shift Model
#' 
#' get Phase I corrected charting constant
#' 
#' @aliases getCC.XBAR
#' @aliases getCC
#' 
#' @param m number of subgroups when the data are subgrouped or number of observations when the data are individual.
#' @param fap0 nominal False Alarm Probabilty in Phase 1
#' @param var.est 'S' - use mean-square-based estimator, 'MR' - use moving-range-based estimator
#' @param ub.cons unbiasing constant 
#' @param method 'exact' - calculate results using the exact method, 'BA' - calculate results using the Bonfferoni approximation 
#' @param interval searching range of charting constants for the exact method 
#' @param nsim number of simulation for the exact method 
#' @param nu degrees of freedom; When var.est = 'S', the degrees of freedom is that of the chi-squared distribution itself for the variance estimator.  When var.est = 'MR', the degrees of freedom is that of the chi-squared distribution approximating to the actual distribution.
#' @param lambda unbiasing constant for the chi-squared distribution approximation. When var.est = 'S', there is no need to do the unbiasing.  When var.est = 'MR', the unbiasing constant needs to be used.
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
                       interval = c(1, 4),
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
                  interval = c(1, 4),
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
