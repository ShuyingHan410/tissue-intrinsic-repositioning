## R3 (revision plan v17): bias and 95% CI coverage of the proposed estimator D-hat (and the per-arm delta-hat), proposed method only.
## Frozen v16 cells stored p-values only, so the estimator is re-run here with the frozen generator (engine.R / engine_v16.R unchanged):
## k = 10 (2/3/3/2 pseudo-mice), M = 64 (plus 27, 331), NB size 1, gate >=30 cells & SE0 < .05, REML + Knapp-Hartung per arm, Welch df contrast.
## Truth: nominal D (theta solved on the arm's cells, u = 0) AND the realised random-effects mean = mean over the arm's slices of
## delta_exp(r_s, theta + u_s) (the estimand of a random-effects model when slice gradients are nudged).
## Run from simulation/v16/:  R_N=500 CORES=5 Rscript ../v17_revision/r3_bias_coverage.R
suppressMessages({source("engine_v16.R"); library(parallel)})
R_N <- as.integer(Sys.getenv("R_N", "500")); CORES <- as.integer(Sys.getenv("CORES", "5")); OUT <- "../v17_revision"
x <- load_substrate("../../.."); pool2 <- substrate_pool(x, "Fibro 2"); d_real <- make_real_arms(x)
pool_ci <- function(tb) { tb <- tb %>% dplyr::filter(is.finite(delta), is.finite(se), se > 0, se < SE_MAX, n >= MIN_CELLS); k <- nrow(tb)
  if (k < 2) return(list(mu = NA, se = NA, lo = NA, hi = NA, k = k))
  fit <- try(rma(yi = tb$delta, sei = tb$se, method = "REML", test = "knha"), silent = TRUE); if (inherits(fit, "try-error")) return(list(mu = NA, se = NA, lo = NA, hi = NA, k = k))
  list(mu = as.numeric(fit$beta), se = fit$se, lo = fit$ci.lb, hi = fit$ci.ub, k = k) }
CELLS <- rbind(
  expand.grid(sub = "pseudo", scen = "gained", D = c(.025, .05, .075), tau = c(0, .05, .15, .30), M = 64, stringsAsFactors = FALSE),
  expand.grid(sub = "pseudo", scen = "null1",  D = 0, tau = c(.05, .15), M = 64, stringsAsFactors = FALSE),
  expand.grid(sub = "pseudo", scen = "gained", D = .05, tau = .05, M = c(27, 331), stringsAsFactors = FALSE),
  expand.grid(sub = "pseudo", scen = "reversed", D = .05, tau = .15, M = 64, stringsAsFactors = FALSE),
  expand.grid(sub = "real",   scen = "gained", D = .05, tau = .15, M = 64, stringsAsFactors = FALSE))
one_rep <- function(i, cl) {
  d <- if (cl$sub == "real") d_real else design_k(pool2, 10)
  rA <- d$r[d$C == "0"]; rB <- d$r[d$C == "1"]
  tc <- switch(cl$scen, null1 = c(0, 0), gained = c(0, solve_theta(if (cl$sub == "real") rB else rA, cl$D)), reversed = c(solve_theta(rA, -cl$D / 2), solve_theta(rB, cl$D / 2)))
  us <- setNames(rnorm(nlevels(d$slice), 0, cl$tau), levels(d$slice))
  dd <- sim_design_u(d, tc, cl$M, 1, u_s = us)
  ## realised estimand: mean over slices of the slice's expected delta at its nudged theta
  sl <- split(seq_len(nrow(d)), d$slice); arm_of <- sapply(sl, function(ix) as.character(d$C[ix[1]]))
  dl <- sapply(names(sl), function(s) delta_exp(d$r[sl[[s]]], tc[as.integer(arm_of[[s]]) + 1] + us[[s]]))
  real0 <- mean(dl[arm_of == "0"]); real9 <- mean(dl[arm_of == "1"])
  tb <- per_slice_table(dd); a0 <- pool_ci(tb %>% dplyr::filter(C == "0")); a9 <- pool_ci(tb %>% dplyr::filter(C == "1"))
  if (!is.finite(a0$se) || !is.finite(a9$se)) return(rep(NA_real_, 14))
  Dh <- a9$mu - a0$mu; se <- sqrt(a0$se^2 + a9$se^2); nu <- se^4 / (a9$se^4 / max(a9$k - 1, 1) + a0$se^4 / max(a0$k - 1, 1)); q <- qt(.975, max(nu, 1))
  c(D_hat = Dh, D_se = se, D_lo = Dh - q * se, D_hi = Dh + q * se, D_nom = cl$D * (cl$scen != "null1"), D_real = real9 - real0,
    d9_hat = a9$mu, d9_lo = a9$lo, d9_hi = a9$hi, d9_nom = tc[2] * 0 + (if (cl$scen == "null1") 0 else if (cl$scen == "reversed") cl$D / 2 else cl$D), d9_real = real9, k0 = a0$k, k9 = a9$k, nu = nu)
}
res <- NULL
for (id in seq_len(nrow(CELLS))) { cl <- CELLS[id, ]; RNGkind("L'Ecuyer-CMRG"); set.seed(20261101 + id); t0 <- Sys.time()
  P <- do.call(rbind, mclapply(seq_len(R_N), one_rep, cl = cl, mc.cores = CORES, mc.set.seed = TRUE)); P <- P[is.finite(P[, 1]), , drop = FALSE]
  s <- data.frame(cl, n = nrow(P),
    bias_nom = mean(P[, "D_hat"] - P[, "D_nom"]), bias_real = mean(P[, "D_hat"] - P[, "D_real"]), sd_D = sd(P[, "D_hat"]), mean_se = mean(P[, "D_se"]),
    cover_nom = mean(P[, "D_lo"] <= P[, "D_nom"] & P[, "D_hi"] >= P[, "D_nom"]), cover_real = mean(P[, "D_lo"] <= P[, "D_real"] & P[, "D_hi"] >= P[, "D_real"]),
    d9_bias_nom = mean(P[, "d9_hat"] - P[, "d9_nom"]), d9_bias_real = mean(P[, "d9_hat"] - P[, "d9_real"]),
    d9_cover_nom = mean(P[, "d9_lo"] <= P[, "d9_nom"] & P[, "d9_hi"] >= P[, "d9_nom"]), d9_cover_real = mean(P[, "d9_lo"] <= P[, "d9_real"] & P[, "d9_hi"] >= P[, "d9_real"]),
    mean_nu = mean(P[, "nu"]), secs = round(as.numeric(Sys.time() - t0, units = "secs")))
  cat(sprintf("DONE %-6s %-8s D=%.3f tau=%.2f M=%3d  n=%d  bias %+.4f (real %+.4f)  SD %.4f  meanSE %.4f  cover %.3f (real %.3f)  d9 cover %.3f (real %.3f)\n",
      cl$sub, cl$scen, cl$D, cl$tau, cl$M, s$n, s$bias_nom, s$bias_real, s$sd_D, s$mean_se, s$cover_nom, s$cover_real, s$d9_cover_nom, s$d9_cover_real)); flush.console()
  res <- rbind(res, s); write.csv(res, file.path(OUT, "r3_bias_coverage.csv"), row.names = FALSE); saveRDS(P, file.path(OUT, sprintf("r3_raw_cell_%02d.rds", id))) }
cat("\nALL DONE\n")
