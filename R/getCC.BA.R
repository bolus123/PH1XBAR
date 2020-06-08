getCC.BA <- function(FAP0, m, nu, ubCons = 1, lambda = 1) {
  ubCons / lambda * qt(((1 - FAP0) ^ (1/m) + 1) / 2, nu)
}
