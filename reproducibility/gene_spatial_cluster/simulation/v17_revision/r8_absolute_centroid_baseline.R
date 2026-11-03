## R8 (added 2026-10-31 after external advice): conceptual baseline = ABSOLUTE transcript centroid T_s = sum(Y r)/sum(Y) per slice
## (no subtraction of the population's own mean position), pooled per arm with the SAME REML + Knapp-Hartung machinery and contrasted
## between arms. Purpose: isolate the effect of the population-relative subtraction. Expected: identical to the proposed method when
## the population does not move (Null 1 / gained on matched pseudo-arms), rejects under whole-population movement (Null 3A/3B tilt,
## real D0 Fibro 2 vs D9 IAF 4 geometry) where the proposed method does not. 200 reps/cell, both statistics on the same data.
## Run from simulation/v16/:  R_N=200 CORES=5 Rscript ../v17_revision/r8_absolute_centroid_baseline.R
suppressMessages({source("engine_v16.R"); source(file.path(PILOT, "null3_resampler.R")); library(parallel)})
R_N <- as.integer(Sys.getenv("R_N", "200")); CORES <- as.integer(Sys.getenv("CORES", "5")); OUT <- "../v17_revision"
x <- load_substrate("../../.."); pool2 <- substrate_pool(x, "Fibro 2"); sm <- readRDS(file.path(PILOT, "null3A_smoke.rds")); TL <- sm$TL; d_real <- make_real_arms(x)
slice_abs <- function(r, Y) { tot <- sum(Y); if (tot <= 0) return(c(delta = NA, se = NA)); Neff <- tot^2 / sum(Y^2); v <- var(r) / Neff
  c(delta = sum(Y * r) / tot, se = if (is.finite(v) && v > 0) sqrt(v) else NA) }
m_absolute <- function(d) { tb <- d %>% dplyr::group_by(slice, mouse, C) %>% dplyr::summarise(n = n(), as.data.frame(t(slice_abs(r, Y))), .groups = "drop")
  a0 <- pool_arm(tb %>% dplyr::filter(C == "0")); a9 <- pool_arm(tb %>% dplyr::filter(C == "1")); unname(d1_test(a0, a9)["p"]) }
substrate_for <- function(sub) switch(sub, pseudo = design_k(pool2, 10), tilt = make_pseudo_arms_null3(pool2, TL$w, 4, c(2, 3, 3, 2)), real = d_real)
CELLS <- rbind(
  data.frame(sub = "pseudo", scen = c("null1", "gained"), D = c(0, .05), tau = .05),
  data.frame(sub = "tilt",   scen = c("null3B", "null3A", "gained"), D = c(0, 0, .05), tau = .15),
  data.frame(sub = "real",   scen = c("null1", "gained"), D = c(0, .05), tau = .15), stringsAsFactors = FALSE)
one_rep <- function(i, cl) { d <- substrate_for(cl$sub); rA <- d$r[d$C == "0"]; rB <- d$r[d$C == "1"]
  tc <- switch(cl$scen, null1 = c(0, 0), null3B = c(0, 0), null3A = c(sm$th0A, sm$th9A), gained = c(0, solve_theta(if (cl$sub == "pseudo") rA else rB, cl$D)))
  dd <- sim_design(d, tc, 64, 1, tau_theta = cl$tau); c(proposed = m_proposed(dd)$pval, absolute = m_absolute(dd), pop_shift = mean(rB) - mean(rA)) }
res <- NULL
for (id in seq_len(nrow(CELLS))) { cl <- CELLS[id, ]; RNGkind("L'Ecuyer-CMRG"); set.seed(20261103 + id); t0 <- Sys.time()
  P <- do.call(rbind, mclapply(seq_len(R_N), one_rep, cl = cl, mc.cores = CORES, mc.set.seed = TRUE))
  row <- data.frame(cl, n_rep = nrow(P), proposed_rej = round(100 * mean(P[, "proposed"] < .05, na.rm = TRUE), 1), absolute_rej = round(100 * mean(P[, "absolute"] < .05, na.rm = TRUE), 1),
                    mean_pop_shift = round(mean(P[, "pop_shift"]), 3), secs = round(as.numeric(Sys.time() - t0, units = "secs")))
  cat(sprintf("DONE %-6s %-6s D=%.3f tau=%.2f  proposed %5.1f%%  absolute %5.1f%%  pop shift %+.3f  (%ds)\n", cl$sub, cl$scen, cl$D, cl$tau, row$proposed_rej, row$absolute_rej, row$mean_pop_shift, row$secs)); flush.console()
  res <- rbind(res, row); write.csv(res, file.path(OUT, "r8_absolute_centroid_baseline.csv"), row.names = FALSE) }
cat("ALL DONE\n")
