PH1XBAR <- function(X,
                    cc = NULL,
                    fap0 = 0.1,
                    var.est = c("S", "MR"),
                    ub.option = TRUE,
                    method = c("exact", "BA"),
                    plot.option = TRUE,
                    interval = c(1, 5),
                    nsim = 10000,
                    verbose = FALSE) {
  var.est <- var.est[1]
  method <- method[1]

  m <- dim(X)[1]
  n <- dim(X)[2]

  X.bar <- rowMeans(X)

  X.bar.bar <- mean(X.bar)

  ub.cons <- 1


  if (var.est == "S") {
    nu <- m - 1
    lambda <- 1

    ub.cons <- ifelse(ub.option == TRUE, c4.f(nu), 1)

    sigma.v <- sqrt(var(X.bar)) / ub.cons
  } else if (var.est == "MR") {
    pars <- pars.root.finding(m - 1, 2, lower = 1e-6)

    nu <- pars[1]
    lambda <- pars[2]

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
      lambda = lambda,
      verbose = verbose
    )
  } else {
    cc <- cc
  }

  LCL <- X.bar.bar - cc * sigma.v
  UCL <- X.bar.bar + cc * sigma.v

  if (plot.option == TRUE) {
    main.text <- paste("Phase I X-bar Chart for fap0 =", fap0)

    plot(c(1, m), c(LCL, UCL), xaxt = "n", xlab = "Subgroup",
        ylab = "Sample Mean", type = "n", main = main.text)

    axis(side = 1, at = 1:m)

    points(1:m, X.bar, type = "o")
    points(c(-1, m + 2), c(LCL, LCL), type = "l", col = "red")
    points(c(-1, m + 2), c(UCL, UCL), type = "l", col = "red")
    points(c(-1, m + 2), c(X.bar.bar, X.bar.bar), type = "l", col = "blue")
    text(round(m * 0.8), UCL, paste("UCL = ", round(UCL, 4)), pos = 1)
    text(round(m * 0.8), LCL, paste("LCL = ", round(LCL, 4)), pos = 3)
  }

  res <- list(CL = X.bar.bar, var.est = sigma.v * ub.cons, ub.cons = ub.cons,
              cc = cc, m = m, nu = nu, lambda = lambda, LCL = LCL, UCL = UCL,
              CS = X.bar)
  invisible(res)
}
