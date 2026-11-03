## Proposed-method-only previews at the v16 reference design (k = 10, M = 64, size 1):
##   Experiment 2 grid: gained/reversed x |D| in {.025,.05,.075} x tau in {0,.05,.15,.30}
##   Experiment 8 slice: k = 20, gained .05, tau in {0,.05,.15,.30}
##   Experiment 4B:      REML+KH vs fixed-effect, null + gained, tau in {0,.18}, k = 10
## Nominal and matched power (matched cutoff = 5th pct of this run's own null p at the same k, tau).
## Output: power_preview.csv, ablB_k10.csv, power_preview.rds
suppressMessages({source("engine_v16.R"); library(parallel)})
R_N <- as.integer(Sys.getenv("R_N", "200")); M_MED <- 64; SIZE <- 1
CORES <- max(1L, min(6L, detectCores() - 1L))
x <- load_substrate(); pool <- substrate_pool(x)
taus <- c(0, .05, .15, .30); Ds <- c(.025, .05, .075)

run_cell <- function(k, tau, pattern, D, fixed = FALSE, seed) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  theta_c <- switch(pattern,
    null     = c(0, 0),
    gained   = c(0, solve_theta(pool$r, D)),
    reversed = c(solve_theta(pool$r, -D / 2), solve_theta(pool$r, D / 2)))
  unlist(mclapply(seq_len(R_N), mc.cores = CORES, function(i) {
    d <- design_k(pool, k); dd <- sim_design(d, theta_c, M_MED, SIZE, tau_theta = tau)
    m_proposed(dd, fixed = fixed)$pval }))
}
## null p-values per (k, tau) for matched cutoffs -- 50 genes x R_N to mirror the A2 rule
null_p <- list()
for (k in c(10, 20)) for (tau in taus) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(30000 + 100 * k + round(100 * tau))
  null_p[[sprintf("%d|%.2f", k, tau)]] <- unlist(mclapply(seq_len(R_N), mc.cores = CORES, function(i) {
    d <- design_k(pool, k); vapply(1:50, function(g) m_proposed(sim_design(d, c(0, 0), M_MED, SIZE, tau_theta = tau))$pval, numeric(1)) }))
  cat(sprintf("null k=%d tau=%.2f  TypeI %.4f  cutoff %.4f\n", k, tau,
              mean(null_p[[sprintf("%d|%.2f", k, tau)]] < .05, na.rm = TRUE), quantile(null_p[[sprintf("%d|%.2f", k, tau)]], .05, na.rm = TRUE)))
}
## Experiment 2 preview (k = 10) and Experiment 8 preview (k = 20, D = .05)
grid <- rbind(expand.grid(k = 10, tau = taus, pattern = c("gained", "reversed"), D = Ds, stringsAsFactors = FALSE),
              expand.grid(k = 20, tau = taus, pattern = "gained", D = .05, stringsAsFactors = FALSE))
res <- list()
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]; p <- run_cell(g$k, g$tau, g$pattern, g$D, seed = 40000 + i)
  cut <- quantile(null_p[[sprintf("%d|%.2f", g$k, g$tau)]], .05, na.rm = TRUE)
  res[[i]] <- data.frame(g, nominal = mean(p < .05, na.rm = TRUE), matched = mean(p < cut, na.rm = TRUE), cutoff = cut, n = sum(is.finite(p)))
  cat(sprintf("k=%2d tau=%.2f %-8s D=%.3f  nominal %.3f  matched %.3f\n", g$k, g$tau, g$pattern, g$D, res[[i]]$nominal, res[[i]]$matched))
}
pw <- bind_rows(res); write_csv(pw, "power_preview.csv")
## Experiment 4B at k = 10: REML+KH vs fixed
ab <- list()
for (tau in c(0, .18)) for (pat in c("null", "gained")) for (D in if (pat == "null") 0 else Ds) for (fx in c(FALSE, TRUE)) {
  p <- run_cell(10, tau, pat, D, fixed = fx, seed = 50000 + round(1000 * tau) + 10 * D * 100 + fx)
  ab[[length(ab) + 1]] <- data.frame(tau = tau, pattern = pat, D = D, pooling = if (fx) "fixed" else "REML+KH", reject = mean(p < .05, na.rm = TRUE))
  cat(sprintf("4B tau=%.2f %-6s D=%.3f %-8s reject %.3f\n", tau, pat, D, if (fx) "fixed" else "REML+KH", mean(p < .05, na.rm = TRUE)))
}
ab <- bind_rows(ab); write_csv(ab, "ablB_k10.csv")
saveRDS(list(null_p = null_p, power = pw, ablB = ab), "power_preview.rds")
cat("DONE 06_power_preview\n")
