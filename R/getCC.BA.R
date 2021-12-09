getCC.ba <- function(fap0, m, nu, ub.cons = 1, lambda = 1) {
  ub.cons / lambda * qt(((1 - fap0)^(1 / m) + 1) / 2, nu)
}
