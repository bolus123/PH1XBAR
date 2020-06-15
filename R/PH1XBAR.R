PH1XBAR <- function(
  X
  ,cc = NULL
  ,FAP0 = 0.1
  ,var.est = c('S', 'MR')
  ,ub.option = TRUE
  ,method = c('exact', 'BA')
  ,plot.option = TRUE
  ,interval = c(1, 5)
  ,nsim = 10000 
  ,seed = 12345
) {
  
  var.est <- var.est[1]
  method <- method[1]
  
  m <- dim(X)[1]
  n <- dim(X)[2]
  
  X.bar <- rowMeans(X)
  
  X.bar.bar <- mean(X.bar)
  
  ubCons <- 1
  
  
  if (var.est == 'S'){
    
    nu <- m - 1
    lambda <- 1
    
    ubCons <- ifelse(ub.option == TRUE, c4.f(nu), 1)
    
    sigma.v <- sqrt(var(X.bar)) / ubCons
    
  } else if (var.est == 'MR') {
    
    pars <- pars.root.finding(m - 1, 2, lower = 1e-6) 
    
    nu <- pars[1]
    lambda <- pars[2]
    
    ubCons <- ifelse(ub.option == TRUE, 1.128, 1)
    
    sigma.v <- mean(abs(diff(X.bar))) / ubCons
    
  }
  
  
  if (is.null(cc)) {
    cc <- getCC(
      FAP0 = FAP0
      ,m = m
      ,var.est = var.est
      ,ubCons = ubCons 
      ,method = method
      ,interval = interval
      ,nsim = nsim
      ,seed = seed
      ,nu = nu 
      ,lambda = lambda
      
    )
  } else {
    
    cc <- cc
    
  }
  
  LCL <- X.bar.bar - cc * sigma.v
  UCL <- X.bar.bar + cc * sigma.v
  
  if (plot.option == TRUE){
    
	main.text <- paste('Phase I X-bar Chart for FAP0 =', FAP0)
	
    plot(c(1, m), c(LCL, UCL), xaxt = "n", xlab = 'Subgroup', ylab = 'Sample Mean', type = 'n', main = main.text)
    
    axis(side = 1, at = 1:m)
    
    points(1:m, X.bar, type = 'o')
    points(c(-1, m + 2), c(LCL, LCL), type = 'l', col = 'red')
    points(c(-1, m + 2), c(UCL, UCL), type = 'l', col = 'red')
    points(c(-1, m + 2), c(X.bar.bar, X.bar.bar), type = 'l', col = 'blue')
    text(round(m * 0.8), UCL, paste('UCL = ', round(UCL, 4)), pos = 1)
    text(round(m * 0.8), LCL, paste('LCL = ', round(LCL, 4)), pos = 3)
    
    
  }
  
  list(CL = X.bar.bar, var.est = sigma.v * ubCons, ubCons = ubCons, cc = cc, m = m, nu = nu, lambda = lambda, LCL = LCL, UCL = UCL, CS = X.bar)
  
}
