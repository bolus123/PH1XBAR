ub.cons.f <- function(nu, ub.option) {

    if (ub.option == 'c4') {
        ub.cons <- c4.f(nu)
    } else if (ub.option == 'd2') {
        ub.cons <- d2.f(2)
    } else {
        ub.cons <- 1
    }
        
    ub.cons
        
}