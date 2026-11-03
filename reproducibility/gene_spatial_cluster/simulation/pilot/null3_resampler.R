## Null 3: moment-balanced population resampler (plan v12 par.1702).
## Weights w_i = exp(g1*r_i + g2*r_i^2) solved jointly so the weighted
## population hits mean = mu_target while keeping SD(r) = SD of the source
## pool:  E_w[r] = mu*,  E_w[r^2] = sigma0^2 + mu*^2.
solve_gamma_tilt <- function(r, mu_target, sd_target = sd(r)) {
  m2_target <- sd_target^2 + mu_target^2
  obj <- function(g) {
    e <- g[1] * r + g[2] * r^2; e <- e - max(e); w <- exp(e)
    m1 <- sum(w * r) / sum(w); m2 <- sum(w * r^2) / sum(w)
    (m1 - mu_target)^2 + (m2 - m2_target)^2
  }
  best <- NULL
  for (st in list(c(0, 0), c(5, -5), c(10, -10), c(-5, 5))) {
    o <- optim(st, obj, method = "BFGS", control = list(maxit = 500, reltol = 1e-14))
    if (is.null(best) || o$value < best$value) best <- o
  }
  g <- best$par
  e <- g[1] * r + g[2] * r^2; e <- e - max(e); w <- exp(e)
  m1 <- sum(w * r) / sum(w); m2 <- sum(w * r^2) / sum(w)
  list(gamma = g, w = w / sum(w),
       ach_mean = m1, ach_sd = sqrt(max(m2 - m1^2, 0)),
       target_mean = mu_target, target_sd = sd_target,
       converged = best$value < 1e-10)
}

## Pseudo-arm substrate where arm 1 cells are drawn from the tilted population.
make_pseudo_arms_null3 <- function(pool, w, mice_per_arm = 4,
                                   slices_per_mouse = c(1, 1, 1, 2)) {
  out <- list(); k <- 0L
  for (arm in c(0, 1)) for (m in seq_len(mice_per_arm)) {
    for (j in seq_len(slices_per_mouse[m])) {
      k <- k + 1L
      n <- sample(pool$sizes, 1)
      rr <- if (arm == 0) sample(pool$r, n, replace = TRUE)
            else sample(pool$r, n, replace = TRUE, prob = w)
      out[[k]] <- data.frame(
        slice = sprintf("A%d_M%d_S%d", arm, m, j),
        mouse = sprintf("A%d_M%d", arm, m),
        C = arm, r = rr)
    }
  }
  d <- dplyr::bind_rows(out)
  d$slice <- factor(d$slice); d$mouse <- factor(d$mouse); d$C <- factor(d$C)
  d
}
