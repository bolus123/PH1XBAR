# moment of W

W.moment <- function(a, n) { # a is order of moment and n is sample size
  
    W.f <- function(w, n) { # n is sample size
  
        integrand <- function(x, w) {
    
            (pnorm(w + x) - pnorm(x)) ^ (n - 2) * dnorm(w + x) * dnorm(x)
    
        }
  
            n * (n - 1) * integrate(integrand, -Inf, Inf, w = w)$value
  
    }
  
    integrand <- function(w, a, n) {
    
        w ^ a * W.f(w, n)
    
    }
  
    integrand <- Vectorize(integrand, 'w')
  
    integrate(integrand, 0, Inf, a = a, n = n)$value
  
}