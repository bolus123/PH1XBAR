c4.f <- function(nu) sqrt(2 / (nu)) / beta((nu) / 2, 1 / 2) * sqrt(pi) # c4.function

#############################################################################################

# moment of W

w.moment <- function(a, n) { # a is order of moment and n is sample size

  W.f <- function(w, n) { # n is sample size

    integrand <- function(x, w) {
      (pnorm(w + x) - pnorm(x))^(n - 2) * dnorm(w + x) * dnorm(x)
    }

    n * (n - 1) * integrate(integrand, -Inf, Inf, w = w)$value
  }

  integrand <- function(w, a, n) {
    w^a * W.f(w, n)
  }

  integrand <- Vectorize(integrand, "w")

  integrate(integrand, 0, Inf, a = a, n = n)$value
}

pars.root.finding <- function(k, n, lower = 1e-6) {
  root.finding <- function(nu, k, n, ratio) {
    ratio / k - nu / 2 / pi * beta(nu / 2, 1 / 2)^2 + 1
  }

  M <- w.moment(1, n)
  M2 <- w.moment(2, n)

  V <- M2 - M^2

  ratio <- V / M^2

  # cat("dn2/Vn = ", 1 / ratio, "\n")

  nu <- uniroot(root.finding, c(lower, k * n), k = k, n = n, ratio = ratio)$root

  lambda <- sqrt(V / k + M^2)

  c(nu, lambda)
}

pars.root.finding <- Vectorize(pars.root.finding, c("k", "n"))

d2.f <- function(n) {
  w.moment(1, n)
}

#############################################################################################
