#' Phase I X-bar control chart with a corrected charting constant
#' 
#' Build a Phase I Shewhart control chart for the variance components model if the data are subgrouped or for the basic Shewhart model if the data are individual. The charting constant is correted by this approach.
#' @param X input and it must be a matrix (m by n) or a vector (m by 1)
#' @param cc nominal Phase I charting constant. If this is given, the function will not recompute the charting constant. 
#' @param fap0 nominal False Alarm Probabilty in Phase 1 
#' @param var.est 'S' - use mean-square-based estimator, 'MR' - use moving-range-based estimator
#' @param ub.option TRUE - the standard deviation estimator corrected by a unbiasing constant. For S, it is c4 and for MR, it is d2. FALSE - no unbiasing constant 
#' @param method 'exact' - calculate results using the exact method, 'BA' - calculate results using the Bonfferoni approximation 
#' @param plot.option - draw a plot for the process; TRUE - Draw a plot for the process, FALSE - Not draw a plot for the process
#' @param interval searching range of charting constants for the exact method 
#' @param nsim number of simulation for the exact method 
#' @param transform type of transformation. When transform = 'none', no transformation is performed. When transform = 'boxcox', the Box-Cox transformation is used. When transform = 'yeojohnson', the Yeo-Johnson transformation is used. 
#' @param lambda parameter used in the Box-Cox or Yeo-Johnson transformation.
#' @param standardize Output standardized charting statistics instead of raw ones. When standardize = TRUE, the standardization is used.  When standardize = FALSE, the standardization is not performed.
#' @param verbose print diagnostic information about fap0 and the charting constant during the simulations for the exact method
#' @return CL Object type double - central line
#' @return var.est Object type double - variance estimate
#' @return ub.cons Object type double - unbiasing constant
#' @return cc Object type double - charting constant
#' @return m Object type integer - number of subgroups when X is a matrix or number of observations when X is a vector
#' @return nu Object type integer - degrees of freedom; When var.est = 'S', the degrees of freedom is that of the chi-squared distribution itself for the variance estimator.  When var.est = 'MR', the degrees of freedom is that of the chi-squared distribution approximating to the actual distribution.
#' @return lambda Object type integer - chi-squared unbiasing constant for the chi-squared distribution approximation
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
                    interval = c(1, 4),
                    nsim = 10000, 
                    transform = "none", 
                    lambda = 1, 
                    standardize=FALSE,
                    verbose = FALSE) {
  var.est <- var.est[1]
  method <- method[1]
  
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  

  m <- dim(X)[1]
  n <- dim(X)[2]

  lambda2 <- NA
  if (transform == "boxcox") {
    lambda2 <- lambda
    X1 <- forecast::BoxCox(X, lambda2)
  } else if (transform == "yeojohnson") {
    lambda2 <- lambda
    X1 <- VGAM::yeo.johnson(X, lambda2)
  } else {
    X1 <- X
  }
  
  X.bar <- rowMeans(X1)

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
  
  
  warning_msg <- NA
  comp_msg <- NULL
  
  if (anyNA(CS) | is.na(X.bar.bar) | is.na(LCL)| is.na(UCL)) {
    if (transform != 'none') {
      warning_msg <- paste(warning_msg, "\nThe back transformation is erroneous using the given lambda.  Please try other lambda's or use standardization.\n", sep = '')
    }
    if (anyNA(CS)) {
      comp_msg <- c(comp_msg, 'One or more CS')
    }
    if (is.na(X.bar.bar)) {
      comp_msg <- c(comp_msg, 'CL')
    }
    if (is.na(LCL)) {
      comp_msg <- c(comp_msg, 'LCL')
    }
    if (is.na(UCL)) {
      comp_msg <- c(comp_msg, 'UCL')
    }
    
    comp_msg <- paste(comp_msg, sep = ', ')
    comp_msg <- paste(comp_msg, " cannot be back-transformed.", sep = '')
    warning_msg <- paste(warning_msg, comp_msg, sep = '')
    warning(warning_msg)
  }
  
  if (plot.option == TRUE) {
    if (n > 1) {
      main.text <- paste("Phase I X-bar Chart")
      plot(c(1, m), c(min(LCL, CS, na.rm = TRUE), max(UCL, CS, na.rm = TRUE)), xaxt = "n", xlab = "Subgroup",
        ylab = "Charting Statistic", type = "n", main = main.text)
    } else {
      main.text <- paste("Phase I Individual Chart")
      plot(c(1, m), c(min(LCL, CS, na.rm = TRUE), max(UCL, CS, na.rm = TRUE)), xaxt = "n", xlab = "Observation",
        ylab = "Charting Statistic", type = "n", main = main.text)
    }
    

    axis(side = 1, at = 1:m)

    points(1:m, CS, type = "o")
    points(c(-1, m + 2), c(LCL, LCL), type = "l", col = "red")
    points(c(-1, m + 2), c(UCL, UCL), type = "l", col = "red")
    if (standardize) {
      points(c(-1, m + 2), c(0, 0), type = "l", col = "blue")
    } else {
      points(c(-1, m + 2), c(X.bar.bar, X.bar.bar), type = "l", col = "blue")
    }
    text(round(m * 0.8), UCL, paste("UCL = ", round(UCL, 4)), pos = 1)
    text(round(m * 0.8), LCL, paste("LCL = ", round(LCL, 4)), pos = 3)
  }

  res <- list(CL = X.bar.bar, var.est = sigma.v * ub.cons, ub.cons = ub.cons,
              cc = cc, m = m, nu = nu, lambda.chi.approx = lambda1, LCL = LCL, UCL = UCL,
              CS = CS, transform = transform, lambda = lambda2, standardize = standardize)
  invisible(res)
}
