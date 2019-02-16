getCC.mvt <- function(
                 m
                 ,nu = m - 1
                 ,FAP = 0.1
                 ,off.diag = -1/(m - 1)
                 ,ub.option = TRUE
                 ,maxiter = 10000

){

    alternative = '2-sided'                                                   #turn off the alternative

    corr.P <- corr.f(m = m, off.diag = off.diag)                              #get correlation matrix with equal correlations

    pu <- 1 - FAP

    ub.cons <- ifelse(ub.option == TRUE, 1, ub.cons.f(nu, 'c4'))

    L <- ifelse(
            alternative == '2-sided',
            qmvt(pu, df = nu, sigma = corr.P, tail = 'both.tails', maxiter = maxiter)$quantile,
            qmvt(pu, df = nu, sigma = corr.P, maxiter = maxiter)$quantile
        )
                                                      #get L by multivariate T



    c.i <- L * ub.cons * sqrt((m - 1) / m)             #get c.i

    list(l = L, c.i = c.i)

}