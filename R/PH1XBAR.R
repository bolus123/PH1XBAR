PH1XBAR <- function(
			X
            ,c.i = NULL
			,FAP = 0.1
			,off.diag = -1/(m - 1)
			#,alternative = '2-sided'
            ,model = 'ANOVA-based'
            ,ub.option = 'c4'
			,plot.option = TRUE
			,maxiter = 10000
			,method = 'direct'
			,indirect.interval = c(1, 7)
			,indirect.subdivisions = 100L
			,indirect.tol = .Machine$double.eps^0.25
) {

    alternative = '2-sided'                                                   #turn off the alternative


	m <- dim(X)[1]
	n <- dim(X)[2]

	X.bar <- rowMeans(X)

	X.bar.bar <- mean(X)

    if (model == 'ANOVA-based') {

        nu <- m - 1

        ub.cons <- ub.cons.f(nu, ub.option)

        sigma.v <- sqrt(var(X.bar)) / ub.cons

    } else if (model == 'basic') {

        nu <- m * (n - 1)

        ub.cons <- ub.cons.f(nu, ub.option)

        sigma.v <- sqrt(sum(apply(X, 1, var)) / m) / ub.cons / sqrt(n)

    } else {

        stop("need to specify whether it is based on the ANOVA-based model or others")

    }

    if (is.null(c.i)) {
        c.i <- PH1.get.cc(
                m = m
                ,nu = nu
                ,FAP = FAP
                ,off.diag = off.diag
                #,alternative = alternative
                ,ub.option = ub.option
                ,maxiter = maxiter
                ,method = method
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
