pars.root.finding <- function(k, n, lower = 1e-6) {
  
    root.finding <- function(nu, k, n, ratio) {
    
        ratio / k - nu / 2 / pi * beta(nu / 2, 1 / 2) ^ 2 + 1
    
    }
  
    M <- W.moment(1, n)
    M2 <- W.moment(2, n)
  
    V <- M2 - M ^ 2
  
    ratio <- V / M ^ 2
  
    cat('dn2/Vn = ', 1 / ratio, '\n')
  
    nu <- uniroot(root.finding, c(lower, k * n), k = k, n = n, ratio = ratio)$root
  
    lambda <- sqrt(V / k + M ^ 2)
  
    c(nu, lambda)
  
}

pars.root.finding <- Vectorize(pars.root.finding, c('k', 'n'))