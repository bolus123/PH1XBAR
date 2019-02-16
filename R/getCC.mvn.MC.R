getCC.mvn.MC <- function(
    m
    ,nu
    ,FAP = 0.1
    #,Phase1 = TRUE
    ,off.diag = -1/(m - 1)
    #,alternative = '2-sided'
    ,var.est = 'MSE'
    ,ub.option = TRUE
    ,ub.lower = 1e-6
    ,sim.X = 10000
    ,sim.Y = 10000
    ,interval = c(1, 7)
    ,maxiter = 10000
    ,tol = 1e-2 
){
    
    require(mvtnorm)
    
    rmvn.MC <- function(sim, sigma) {
  
        rmvnorm(sim, sigma = sigma)
  
    }

    jointPDF.mvn.chisq.MC <- function(
        Y
        ,c.i
        ,m
        ,nu
        ,sigma
        ,lambda = 1
        #,alternative = '2-sided'
        ,ub.cons = 1
        ,X
    ){
  
        alternative = '2-sided'                                                   #turn off the alternative
  
        s <- length(Y)
  
        L <- c.i / sqrt((m - 1) / m * nu) * sqrt(Y) / ub.cons * lambda
  
        dpp <- lapply(
            1:s,
            function(i){
      
                n <- dim(X)[1]
      
                ifelse(
                    alternative == '2-sided',
                    mean(rowSums(-L[i] < X & X < L[i]) == m),
                    mean(rowSums(-Inf < X & X < L[i]) == m)
                )
      
            }
    
        )
  
        dpp <- unlist(dpp)
  
        dpp 

    }


    Root.mvn.f.MC <- function(
        c.i
        , m
        , nu
        , sigma
        , lambda = 1
        , pu
        #, alternative = '2-sided'
        , ub.cons = 1
        , X
        , Y
    ){
        alternative = '2-sided'                                                   #turn off the alternative

        pp <- mean(
            jointPDF.mvn.chisq.MC(
                Y = Y,
                c.i = c.i,
                m = m,
                nu = nu,
                sigma = sigma,
                lambda = lambda,
                #alternative = alternative,
                ub.cons = ub.cons,
                X = X
            )
        )
  
        cat('ci:', c.i, ' and diff:', pu - pp, '\n')
  
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
  
    if (var.est == 'MSE') {
    
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
  
    X <- rmvn.MC(sim = sim.X, sigma = corr.P)
    Y <- rchisq(sim.Y, nu)
  
    pu <- 1 - FAP
  
    c.i <- uniroot(
        Root.mvn.f.MC,
        interval = interval,
        m = m,
        nu = nu,
        #Y = Y,
        sigma = corr.P,
        lambda = lambda,
        pu = pu,
        #alternative = alternative,
        ub.cons = ub.cons,
        X = X,
        Y = Y,
        tol = tol,
        maxiter = maxiter
    )$root
  
    L <- c.i / ub.cons * sqrt(m / (m - 1)) * lambda
  
    list(l = L, c.i = c.i)
  
}