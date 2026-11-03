## =====================================================================
## v16 engine extensions. Sources pilot/engine.R unchanged and ADDS:
##   sim_slice_h / delta_exp_h   arbitrary weight vector h (Exp 11 shapes)
##   h_local, h_bimodal, solve_local, solve_bimodal
##   sim_design_u                sim_design with caller-supplied slice nudges
##                               u_s and mouse effects b_m (Exp 6b, Exp 9 B5)
##   draw_block_nudges           block-correlated u_gs for the FDR mixture
##   draw_mouse_effects          normal / t3 / one-outlier mouse effects
## No function in engine.R is redefined.
## =====================================================================
## Run every v16 script from the v16/ directory (Rscript <script>.R); pilot/ is a sibling.
PILOT <- Sys.getenv("PILOT", "../pilot")
suppressMessages(source(file.path(PILOT, "engine.R")))
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- arbitrary-shape generator ----------------------------------------
sim_slice_h <- function(h, M, size = 2) {
  mu <- M * h / sum(h)
  rnbinom(length(h), mu = mu, size = size)
}
delta_exp_h <- function(r, h) sum(h * r) / sum(h) - mean(r)

## Localised shift: a fraction pi of expected mass moves into [lo, hi]; the rest is flat.
## h sums to 1 by construction. delta(pi) = pi * (mean(r | window) - mean(r)), closed form.
h_local <- function(r, pi, lo = 0.75, hi = 0.85) {
  w <- as.numeric(r >= lo & r <= hi); n <- length(r); nw <- sum(w)
  stopifnot(nw > 0)
  (1 - pi) / n + pi * w / nw
}
solve_local <- function(r, target, lo = 0.75, hi = 0.85) {
  w <- r >= lo & r <= hi
  target / (mean(r[w]) - mean(r))
}

## Bimodal: h = 1 + (kappa/3) phi(r; .25, .05) + kappa phi(r; .75, .05). One parameter.
h_bimodal <- function(r, kappa, m1 = 0.25, m2 = 0.75, s = 0.05, ratio = 1/3) {
  1 + kappa * ratio * dnorm(r, m1, s) / dnorm(0, 0, s) +
      kappa *         dnorm(r, m2, s) / dnorm(0, 0, s)
}
solve_bimodal <- function(r, target, lo = 0, hi = 50, ...) {
  if (abs(target) < 1e-12) return(0)
  f <- function(k) delta_exp_h(r, h_bimodal(r, k, ...)) - target
  if (f(lo) * f(hi) > 0) return(NA_real_)
  uniroot(f, c(lo, hi), tol = 1e-10)$root
}

## ---- design simulation with caller-supplied random effects ------------
## u_s : named numeric, one per slice level (slice nudge on theta scale)
## b_m : named numeric, one per mouse level
## shape: "exp" (default, theta) | "local" | "bimodal"; for the latter two,
##        theta_c is interpreted as the shape parameter (pi or kappa) per arm.
sim_design_u <- function(d, theta_c, M, size = 2, u_s = NULL, b_m = NULL,
                         shape = c("exp", "local", "bimodal")) {
  shape <- match.arg(shape)
  Mv <- rep_len(M, 2)
  sls <- levels(d$slice); mice <- levels(d$mouse)
  if (is.null(u_s)) u_s <- setNames(rep(0, length(sls)), sls)
  if (is.null(b_m)) b_m <- setNames(rep(0, length(mice)), mice)
  d$Y <- NA_integer_
  for (s in sls) {
    i <- which(d$slice == s)
    arm <- as.integer(as.character(d$C[i[1]])) + 1L
    m <- as.character(d$mouse[i[1]])
    if (shape == "exp") {
      th <- theta_c[arm] + b_m[[m]] + u_s[[s]]
      d$Y[i] <- sim_slice(d$r[i], th, Mv[arm], size)
    } else {
      ## shape parameter is not nudged on the theta scale; heterogeneity for
      ## non-exp shapes is applied multiplicatively to the shape parameter
      par <- theta_c[arm] * exp(b_m[[m]] + u_s[[s]])
      h <- if (shape == "local") h_local(d$r[i], par) else h_bimodal(d$r[i], par)
      d$Y[i] <- sim_slice_h(h, Mv[arm], size)
    }
  }
  d
}

## ---- block-correlated slice nudges for the FDR mixture ----------------
## returns G x S matrix of u_gs with marginal sd tau, within-block correlation rho
draw_block_nudges <- function(G, slices, tau, rho = 0.7, block = 20) {
  S <- length(slices)
  nb <- ceiling(G / block)
  blk <- rep(seq_len(nb), each = block)[seq_len(G)]
  W <- matrix(rnorm(nb * S, 0, tau), nb, S)            # shared per block
  E <- matrix(rnorm(G * S, 0, tau), G, S)              # gene-specific
  U <- rho * W[blk, , drop = FALSE] + sqrt(1 - rho^2) * E
  dimnames(U) <- list(paste0("g", seq_len(G)), slices)
  U
}

## ---- mouse-effect distributions (Exp 9, B5) ----------------------------
draw_mouse_effects <- function(mice, sigma, kind = c("normal", "t3", "outlier")) {
  kind <- match.arg(kind)
  n <- length(mice)
  b <- switch(kind,
    normal  = rnorm(n, 0, sigma),
    t3      = sigma * rt(n, df = 3) / sqrt(3),         # variance-matched
    outlier = { v <- rep(0, n); v[sample(n, 1)] <- 3 * sigma; v })
  setNames(b, mice)
}

## ---- convenience: pseudo-arm designs by k ------------------------------
design_k <- function(pool, k = c(5, 10, 20)) {
  k <- match.arg(as.character(k), c("5", "10", "20"))
  switch(k,
    "5"  = make_pseudo_arms(pool, mice_per_arm = 4, slices_per_mouse = c(1, 1, 1, 2)),
    "10" = make_pseudo_arms(pool, mice_per_arm = 4, slices_per_mouse = c(2, 3, 3, 2)),
    "20" = make_pseudo_arms(pool, mice_per_arm = 8, slices_per_mouse = rep(c(2, 3, 3, 2), 2)))
}
