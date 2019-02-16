corr.f <- function(m, off.diag = - 1 / (m - 1)){
                                                                                   #correlation matrix
    crr <- diag(m)
    crr[which(crr == 0)] <- off.diag

    crr

}