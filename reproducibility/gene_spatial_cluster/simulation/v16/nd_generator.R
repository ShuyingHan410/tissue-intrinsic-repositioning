## Non-directional spatial change: h9(r) = 1 + kappa*(r - c)^2 with c solved
## PER SLICE so the expected population-relative delta is exactly zero.
## f(c) = sum(h r)/sum(h) - mean(r) is > 0 at c = min(r) (h increasing in r)
## and < 0 at c = max(r) (h decreasing), so a root always exists in range(r).
solve_c_nd <- function(r, kappa) {
  f <- function(cc) { h <- 1 + kappa * (r - cc)^2; sum(h * r) / sum(h) - mean(r) }
  uniroot(f, range(r), tol = 1e-14)$root
}
delta_exp_nd <- function(r, kappa, cc) {
  h <- 1 + kappa * (r - cc)^2
  sum(h * r) / sum(h) - mean(r)
}
sim_slice_nd <- function(r, kappa, M, size = 2) {
  cc <- solve_c_nd(r, kappa)
  h <- 1 + kappa * (r - cc)^2
  mu <- M * h / sum(h)
  rnbinom(length(r), mu = mu, size = size)
}
## achieved profile change of h vs flat, as domain-share difference (K domains)
nd_domain_diff <- function(r, kappa, cc, K = 5) {
  br <- seq(0, 1, length.out = K + 1)
  dom <- cut(r, br, include.lowest = TRUE, labels = FALSE)
  h <- 1 + kappa * (r - cc)^2
  p1 <- tapply(h, factor(dom, seq_len(K)), sum); p1[is.na(p1)] <- 0
  p0 <- tabulate(dom, K)
  max(abs(p1 / sum(p1) - p0 / sum(p0)))
}
