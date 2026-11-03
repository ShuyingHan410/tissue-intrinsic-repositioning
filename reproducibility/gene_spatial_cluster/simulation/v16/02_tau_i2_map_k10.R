## v16 §0a A5: tau_theta -> realised median I2 at k = 10 (null genes, proposed method).
## Mirrors pilot/fix_ablB_tau.R but at the v13 primary design. Output: tau_i2_k10.csv
suppressMessages({source("engine_v16.R"); library(parallel)})
R_N <- as.integer(Sys.getenv("R_N", "100")); M_MED <- 64
CORES <- max(1L, min(6L, detectCores() - 1L))
x <- load_substrate(); pool <- substrate_pool(x)
taus <- c(0, .05, .10, .15, .20, .25, .30)
out <- list()
for (k in c(10, 5)) for (tau in taus) {
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260915 + 1000 * k + round(100 * tau))
  v <- unlist(mclapply(seq_len(R_N), mc.cores = CORES, function(i) {
    d <- design_k(pool, k)
    dd <- sim_design(d, c(0, 0), M_MED, 2, tau_theta = tau)
    r <- m_proposed(dd); c(r$I2_0, r$I2_9)
  }))
  out[[length(out) + 1]] <- data.frame(k = k, tau = tau, I2_median = median(v, na.rm = TRUE),
                                      I2_p25 = quantile(v, .25, na.rm = TRUE), I2_p75 = quantile(v, .75, na.rm = TRUE))
  cat(sprintf("k=%2d tau=%.2f  realised I2 median %.1f  IQR [%.1f, %.1f]\n", k, tau,
              median(v, na.rm = TRUE), quantile(v, .25, na.rm = TRUE), quantile(v, .75, na.rm = TRUE)))
}
tab <- do.call(rbind, out); write_csv(tab, "tau_i2_k10.csv")
cat("DONE 02_tau_i2_map_k10\n")
