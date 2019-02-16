getCC <- function(
            m
            ,nu = m - 1
            ,FAP = 0.1
            ,off.diag = -1/(m - 1)
            #,alternative = '2-sided'
            ,var.est = 'MS'
            ,ub.option = TRUE
            ,maxiter = 10000
            ,ub.lower = 1e-6
            ,indirect.interval = c(1, 7)
            ,indirect.subdivisions = 100L
            ,indirect.tol = .Machine$double.eps^0.25


){
    alternative = '2-sided'                                                   #turn off the alternative
                                                  #The purpose of this function is to obtain L and K
                                                    #by multivariate T or multivariate normal.
                                                    #Multivariate normal is so time-consuming
                                                    #that I do not recommend.


    is.int <- ifelse(nu == round(nu), 1, 0)
    
    
    if (is.int == 1) {
    
        cat('Nu is an integer. The direct method is being employed.', '\n')
        
        getCC.mvt(
                    m = m
                    ,nu = nu
                    ,FAP = FAP
                    #,Phase1 = Phase1
                    ,off.diag = off.diag
                    #,alternative = alternative
                    ,ub.option = ub.option
                    ,maxiter = maxiter
        )
    
    } else {
    
        cat('Nu is a non-integer. The indirect method is being employed.', '\n')
    
        getCC.mvn(
            m = m
            ,nu = nu
            ,FAP = FAP
            #,Phase1 = Phase1
            ,off.diag = off.diag
            #,alternative = alternative
            ,var.est = var.est
            ,ub.option = ub.option
            ,ub.lower = ub.lower
            ,interval = indirect.interval
            #,maxsim = indirect.maxsim
            ,subdivisions = indirect.subdivisions
            ,maxiter = maxiter
            ,tol = indirect.tol
        )
        
    }


}