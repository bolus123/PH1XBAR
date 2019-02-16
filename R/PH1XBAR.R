PH1XBAR <- function(
    X
    ,c.i = NULL
		,FAP = 0.1
		,off.diag = -1/(m - 1)
    ,var.est = 'MS'
    ,ub.option = TRUE
		,plot.option = TRUE
		,maxiter = 10000
    ,ub.lower = 1e-6
		,indirect.interval = c(1, 7)
		,indirect.subdivisions = 100L
		,indirect.tol = .Machine$double.eps^0.25
) {

    alternative = '2-sided'                                                   #turn off the alternative


	m <- dim(X)[1]
	n <- dim(X)[2]

	X.bar <- rowMeans(X)

	X.bar.bar <- mean(X.bar)

    ub.cons <- 1
    
    if (ub.option == TRUE) {
    
        if (var.est == 'MS'){
        
            nu <- m - 1
        
            ub.cons <- ub.cons.f(nu, 'c4')
            
            sigma.v <- sqrt(var(X.bar)) / ub.cons
        
        } else if (var.est == 'MR') {
        
            pars <- pars.root.finding(m - 1, 2, lower = ub.lower) 
            
            nu <- pars[1]
        
            ub.cons <- ub.cons.f(nu, 'd2')
            
            sigma.v <- mean(abs(X.bar[2:m] - X.bar[1:(m - 1)])) / ub.cons
        
        }
    
    }

    if (is.null(c.i)) {
         c.i <- getCC(
            m = m
            ,nu = nu
            ,FAP = FAP
            ,off.diag = off.diag
            #,alternative = '2-sided'
            ,ub.option = ub.option
            ,maxiter = maxiter
            ,var.est = var.est
            ,ub.lower = ub.lower
            ,indirect.interval = indirect.interval
            ,indirect.subdivisions = indirect.subdivisions
            ,indirect.tol = indirect.tol
        )$c.i
    } else {

            c.i <- c.i

    }

	LCL <- X.bar.bar - c.i * sigma.v
	UCL <- X.bar.bar + c.i * sigma.v

	if (plot.option == TRUE){

		plot(c(1, m), c(LCL, UCL), xaxt = "n", xlab = 'Subgroup', ylab = 'Sample Mean', type = 'n')

		axis(side = 1, at = 1:m)

		points(1:m, X.bar, type = 'o')
		points(c(-1, m + 2), c(LCL, LCL), type = 'l', col = 'red')
		points(c(-1, m + 2), c(UCL, UCL), type = 'l', col = 'red')
		points(c(-1, m + 2), c(X.bar.bar, X.bar.bar), type = 'l', col = 'blue')
		text(round(m * 0.8), UCL, paste('UCL = ', round(UCL, 4)), pos = 1)
		text(round(m * 0.8), LCL, paste('LCL = ', round(LCL, 4)), pos = 3)


	}

	list(CL = X.bar.bar, sigma = sigma.v, PH1.cc = c.i, m = m, nu = nu, LCL = LCL, UCL = UCL, CS = X.bar)

}
