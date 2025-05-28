#' Phase I X-bar control chart with a corrected charting constant
#' 
#' Build a Phase I Shewhart control chart for the variance components model if the data are subgrouped or for the basic Shewhart model if the data are individual. The charting constant is correted by this approach.
#' @param X input and it must be a matrix 
#' @param cc nominal Phase I charting constant. If this is given, the function will not recompute the charting constant. 
#' @param fap0 nominal false Alarm Probabilty in Phase 1 
#' @param var.est 'S' - use mean-square-based estimator, 'MR' - use moving-range-based estimator
#' @param ub.option TRUE - the standard deviation estimator corrected by a unbiasing constant. For S, it is c4 and for MR, it is d2. FALSE - no unbiasing constant 
#' @param method 'exact' - calculate results using the exact method, 'BA' - calculate results using the Bonfferoni approximation 
#' @param plot.option - draw a plot for the process; FALSE - Not draw a plot for the process
#' @param interval searching range of charting constants for the exact method 
#' @param nsim number of simulation for the exact method 
#' @param transform type of transformation. When transform = 'none', no transformation is performed. When transform = 'boxcox', the Box-Cox transformation is used. When transform = 'yeojohnson', the Yeo-Johnson transformation is used. 
#' @param lambda parameter used in the transformation.
#' @param standardize Output standardized data instead of raw data.
#' @param verbose print diagnostic information about fap0 and the charting constant during the simulations for the exact method
#' @return CL Object type double - central line
#' @return var.est Object type double - variance estimate
#' @return ub.cons Object type double - unbiasing constant
#' @return cc Object type double - charting constant
#' @return m Object type integer - number of observations
#' @return nu Object type integer - degrees of freedom
#' @return lambda Object type integer - chi-squared unbiasing constant
#' @return LCL Object type double - lower charting limit
#' @return UCL Object type double - upper charting limit
#' @return CS Object type double - charting statistic
#' @references Champ, C.W., and Jones, L.A. (2004) Designing Phase I X-bar charts with small sample sizes. Quality and Reliability Engineering International. 20(5), 497-510
#' @references Yao, Y., Hilton, C.W., and Chakraborti, S. (2017) Designing Phase I Shewhart X-bar charts: Extended tables and software. Quality and Reliability Engineering International. 33(8), 2667-2672.
#' @references Yao, Y., and Chakraborti, S. (2021). Phase I monitoring of individual normal data: Design and implementation. Quality Engineering, 33(3), 443-456.
#' @references Yao, Y., and Chakraborti, S. (2021). Phase I process monitoring: The case of the balanced one-way random effects model. Quality and Reliability Engineering International, 37(3), 1244-1265.
#' 
#' 
#' @export
#' @examples
#' 
#' 
#' set.seed(12345)
#' 
#' # load the data in the package as an example
#' data(grinder_data)
#' 
#' # An example using a false alarm probability of 0.05, and 10 simulations
#' PH1XBAR(grinder_data, fap0 = 0.05, nsim=10, verbose=TRUE)
#' 
#' 
PH1XBAR <- function(X,
                    cc = NULL,
                    fap0 = 0.05,
                    var.est = c("S", "MR"),
                    ub.option = TRUE,
                    method = c("exact", "BA"),
                    plot.option = TRUE,
                    interval = c(1, 5),
                    nsim = 10000, 
                    transform = "none", 
                    lambda = 1, 
                    standardize=FALSE,
                    verbose = FALSE) {
  var.est <- var.est[1]
  method <- method[1]

  m <- dim(X)[1]
  n <- dim(X)[2]

  lambda2 <- NA
  if (transform == "boxcox") {
    lambda2 <- lambda
    X <- forecast::BoxCox(X, lambda2)
  } else if (transform == "yeojohnson") {
    lambda2 <- lambda
    X <- VGAM::yeo.johnson(X, lambda2)
  }
  
  X.bar <- rowMeans(X)

  X.bar.bar <- mean(X.bar)

  ub.cons <- 1

  if (var.est == "S") {
    nu <- m - 1
    lambda1 <- 1

    ub.cons <- ifelse(ub.option == TRUE, c4.f(nu), 1)

    sigma.v <- sqrt(var(X.bar)) / ub.cons
  } else if (var.est == "MR") {
    pars <- pars.root.finding(m - 1, 2, lower = 1e-6)

    nu <- pars[1]
    lambda1 <- pars[2]

    ub.cons <- ifelse(ub.option == TRUE, 1.128, 1)

    sigma.v <- mean(abs(diff(X.bar))) / ub.cons
  }

  if (is.null(cc)) {
    cc <- getCC.XBAR(
      fap0 = fap0,
      m = m,
      var.est = var.est,
      ub.cons = ub.cons,
      method = method,
      interval = interval,
      nsim = nsim,
      nu = nu,
      lambda = lambda1,
      verbose = verbose
    )
  } 

  if (standardize) {
    CS <- (X.bar - X.bar.bar) / sigma.v
    LCL <- - cc
    UCL <- cc
  } else {
    CS <- X.bar
    LCL <- X.bar.bar - cc * sigma.v
    UCL <- X.bar.bar + cc * sigma.v
    
    if (transform == "boxcox") {
      CS <- forecast::InvBoxCox(CS, lambda2)
      X.bar.bar <- forecast::InvBoxCox(X.bar.bar, lambda2)
      LCL <- forecast::InvBoxCox(LCL, lambda2)
      UCL <- forecast::InvBoxCox(UCL, lambda2)
    } else if (transform == "yeojohnson") {
      CS <- VGAM::yeo.johnson(CS, lambda2, inverse = TRUE)
      X.bar.bar <- VGAM::yeo.johnson(X.bar.bar, lambda2, inverse = TRUE)
      LCL <- VGAM::yeo.johnson(LCL, lambda2, inverse = TRUE)
      UCL <- VGAM::yeo.johnson(UCL, lambda2, inverse = TRUE)
    }
  }
  
  if (plot.option == TRUE) {
    main.text <- paste("Phase I X-bar Chart for fap0 =", fap0)

    plot(c(1, m), c(LCL, UCL), xaxt = "n", xlab = "Subgroup",
        ylab = "Charting Statistic", type = "n", main = main.text)

    axis(side = 1, at = 1:m)

    points(1:m, X.bar, type = "o")
    points(c(-1, m + 2), c(LCL, LCL), type = "l", col = "red")
    points(c(-1, m + 2), c(UCL, UCL), type = "l", col = "red")
    points(c(-1, m + 2), c(X.bar.bar, X.bar.bar), type = "l", col = "blue")
    text(round(m * 0.8), UCL, paste("UCL = ", round(UCL, 4)), pos = 1)
    text(round(m * 0.8), LCL, paste("LCL = ", round(LCL, 4)), pos = 3)
  }

  res <- list(CL = X.bar.bar, var.est = sigma.v * ub.cons, ub.cons = ub.cons,
              cc = cc, m = m, nu = nu, lambda1 = lambda1, LCL = LCL, UCL = UCL,
              CS = CS, transform = transform, lambda = lambda2, standardize = standardize)
  invisible(res)
}
