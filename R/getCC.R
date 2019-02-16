getCC <- function(
            m
            ,nu = m - 1
            ,FAP = 0.1
            ,off.diag = -1/(m - 1)
            #,alternative = '2-sided'
            ,ub.option = TRUE
            ,maxiter = 10000
            ,method = 'direct'
            ,var.est = 'MSE'
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
    
    
    if (method == 'direct') {
    
        if (var.est == 'MSE') {
            if (is.int == 1) {
            
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
                
            } else if (is.int == 0) {
                
                stop('Nu is not an integer. Please use the indirect method instead.')
                
            }
        } else {
        
            stop('The variance estimation must be MSE for the direct method. Please use the indirect method instead.')
        
        }
    
    } else if (method == 'indirect') {
    
        if (is.int == 1 & var.est == 'MSE') {
        
            cat('Nu is an integer. Using the indirect method may slow the computation process down.', '\n')
        } 
              
        
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
    
    } else {
    
        stop('Unknown method. Please select the direct method or the indirect method.')
    
    }


}