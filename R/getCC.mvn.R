getCC.mvn <- function(
                 m
                 ,nu = m - 1
                 ,FAP = 0.1
                 #,Phase1 = TRUE
                 ,off.diag = -1/(m - 1)
                 #,alternative = '2-sided'
                 ,var.est = 'MS'
                 ,ub.option = TRUE
                 ,maxiter = 10000
                 ,ub.lower = 1e-6
                 ,interval = c(1, 7)
                 ,subdivisions = 2000
                 ,tol = 1e-2

){

    jointPDF.mvn.chisq <- function(
                        Y
                        ,c.i
                        ,m
                        ,nu
                        ,sigma
                        ,lambda = 1
                        #,alternative = '2-sided'
                        ,ub.cons = 1)
    {

        alternative = '2-sided'                                                   #turn off the alternative

        s <- length(Y)

        L <- c.i / sqrt((m - 1) / m * nu) * sqrt(Y) / ub.cons * lambda

        dpp <- lapply(
                1:s,
                function(i){

                    LL <- rep(L[i], m)

                    ifelse(
                        alternative == '2-sided',
                        pmvnorm(lower = -LL, upper = LL, sigma = sigma),
                        pmvnorm(lower = rep(-Inf, m), upper = LL, sigma = sigma)
                    )

                }

        )

        dpp <- unlist(dpp)

        dpp * dchisq(Y, df = nu)


    }

    Root.mvn.f <- function(
                        c.i
                        , m
                        , nu
                        , sigma
                        , lambda = 1
                        , pu
                        #, alternative = '2-sided'
                        , ub.cons = 1
                        , subdivisions = 2000
                        , rel.tol = 1e-2)
    {
        alternative = '2-sided'                                                   #turn off the alternative

        pp <- integrate(
                jointPDF.mvn.chisq,
                lower = 0,
                upper = Inf,
                c.i = c.i,
                m = m,
                nu = nu,
                sigma = sigma,
                lambda = lambda,
                #alternative = alternative,
                ub.cons = ub.cons,
                subdivisions = subdivisions,
                rel.tol = rel.tol
            )$value

        pu - pp


    }

    alternative = '2-sided'                                                   #turn off the alternative
                                                         #The purpose of this function is
                                                            #to obtain L and K based on
    #MCMC <- FALSE                                           #multivariate normal.
                                                            #MCMC part is not available now.
    #if (is.null(off.diag)) off.diag <- ifelse(Phase1 == TRUE, - 1 /(m - 1), 1 / (m + 1))

    ub.cons <- 1
    lambda <- 1
    
    if (var.est == 'MS') {
    
        if (ub.option == TRUE) {
        
            ub.cons <- ub.cons.f(nu, 'c4')
            
        }
    
    } else if (var.est == 'MR') {
    
        if (ub.option == TRUE) {
        
            ub.cons <- ub.cons.f(nu, 'd2')
        
        }
        
            nu.lambda <- pars.root.finding(m - 1, 2, lower = ub.lower)
            
            nu <- nu.lambda[1]
            lambda <- nu.lambda[2]
        
        
    
    } else {
    
        stop('The variance estimation is unknown.')
    
    }
    
    cat('nu:', nu, ', lambda:', lambda, '\n')

    corr.P <- corr.f(m = m, off.diag = off.diag)

    pu <- 1 - FAP

    c.i <- uniroot(
            Root.mvn.f,
            interval = interval,
            m = m,
            nu = nu,
            #Y = Y,
            sigma = corr.P,
            lambda = lambda,
            pu = pu,
            #alternative = alternative,
            ub.cons = ub.cons,
            subdivisions = subdivisions,
            tol = tol,
            rel.tol = tol,
            maxiter = maxiter
    )$root
    
    

    L <- c.i / ub.cons * sqrt(m / (m - 1)) * lambda

    list(l = L, c.i = c.i)


}