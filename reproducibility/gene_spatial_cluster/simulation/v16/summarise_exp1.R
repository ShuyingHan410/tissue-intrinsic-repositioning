## Inspect Experiment 1 / 1b (plan v16 §0.5 steps 1-3). Run after all 12 cells exist in exp1_out/.
## Produces:
##   exp1_typeI.csv      Type I per method x cell, cluster-robust MC SE, n_eff, pass flag [3%, 7%]
##   exp1_thresholds.csv c_{m,k,tau} = 5th percentile of each method's own Null 1 p-values (design-specific)
##   exp1_F0.rds         empirical null CDFs F0_{m,k,tau} (sorted Null 1 p-values) for the FDR repair
##   exp1_null2.csv      Null 2 rejection rates -- the "is there a gradient?" vs "did it change?" check
##   exp1_qq.png         QQ of -log10 p against uniform, per method, per cell
suppressMessages({source("engine_v16.R"); library(ggplot2)})
fs <- list.files("exp1_out", "^cell_\\d+[.]rds$", full.names = TRUE)
stopifnot(length(fs) > 0)
cells <- lapply(fs, readRDS)
meth <- c("proposed", "gam_naive", "gam_fs", "despace2", "pbulk_de")

rows <- list(); thr <- list(); F0 <- list(); qq <- list()
for (z in cells) {
  A <- z$P                                       # G0 x 5 x R
  for (m in meth) {
    Pm <- A[, m, ]                               # G0 x R
    if (all(is.na(Pm))) next
    mc <- mc_rate(lapply(seq_len(ncol(Pm)), function(r) Pm[, r] < 0.05))
    rows[[length(rows) + 1]] <- data.frame(z$cell, method = m, t(mc),
      pass_3_7 = mc["alpha"] >= 0.03 & mc["alpha"] <= 0.07, median_p = median(Pm, na.rm = TRUE))
    if (z$cell$null == "null1") {
      p <- sort(as.vector(Pm[is.finite(Pm)]))
      key <- sprintf("%s|k%d|tau%.2f", m, z$cell$k, z$cell$tau)
      thr[[key]] <- data.frame(method = m, k = z$cell$k, tau = z$cell$tau, n_null = length(p),
                               c_0.05 = quantile(p, .05), c_0.01 = quantile(p, .01), c_0.001 = quantile(p, .001))
      F0[[key]] <- p                             # F0(x) = mean(p <= x); stored sorted for ecdf()
      qq[[length(qq) + 1]] <- data.frame(method = m, k = z$cell$k, tau = z$cell$tau,
        expected = -log10(ppoints(length(p))), observed = -log10(p))
    }
  }
}
t1 <- bind_rows(rows); write_csv(t1, "exp1_typeI.csv")
th <- bind_rows(thr);  write_csv(th, "exp1_thresholds.csv")
saveRDS(F0, "exp1_F0.rds")
options(width = 220)
cat("\n## Type I error at nominal 0.05 (Null 1), per method x design\n")
print(as.data.frame(t1 %>% filter(null == "null1") %>% mutate(alpha = sprintf("%.3f [%.3f,%.3f]", alpha, lo, hi)) %>%
  select(exp, k, tau, method, alpha, n_eff, pass_3_7) %>% pivot_wider(names_from = method, values_from = c(alpha, pass_3_7))))
cat("\n## Null 2 rejection rate (positioned identically in both arms) -- large values mean the method tests\n##    'is there a gradient?' rather than 'did it change?' (LOCKED_METHODS.md §3 trap)\n")
n2 <- t1 %>% filter(null == "null2") %>% select(k, tau, method, alpha, lo, hi, n_eff); write_csv(n2, "exp1_null2.csv")
print(as.data.frame(n2 %>% mutate(alpha = sprintf("%.3f [%.3f,%.3f]", alpha, lo, hi)) %>% select(k, tau, method, alpha) %>%
  pivot_wider(names_from = method, values_from = alpha)))
cat("\n## design-specific matched thresholds c_{m,k,tau} (use ONLY at the same k and tau)\n")
print(as.data.frame(th %>% mutate(across(starts_with("c_"), ~signif(.x, 3)))))
## verdicts
cat("\n## verdict per comparator x design\n")
v <- t1 %>% filter(null == "null1", method %in% c("gam_fs", "despace2")) %>%
  mutate(verdict = ifelse(pass_3_7, "nominal: use p < 0.05", sprintf("NON-NOMINAL: use c = %.4f", NA))) %>% select(exp, k, tau, method, alpha, pass_3_7)
v <- v %>% left_join(th %>% select(method, k, tau, c_0.05), by = c("method", "k", "tau")) %>%
  mutate(rule = ifelse(pass_3_7, "nominal p < 0.05", sprintf("matched: p < %.4f", c_0.05)))
print(as.data.frame(v %>% select(exp, k, tau, method, alpha, rule)))
## QQ figure
g <- ggplot(bind_rows(qq) %>% mutate(design = sprintf("k=%d  tau=%.2f", k, tau)), aes(expected, observed, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") + geom_line(linewidth = .6) +
  facet_wrap(~design, ncol = 4) + coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 6)) +
  labs(x = "expected -log10 p (uniform)", y = "observed -log10 p", title = "Experiment 1 / 1b: Null 1 p-value QQ per method") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave("exp1_qq.png", g, width = 11, height = 6.5, dpi = 150)
cat("\nwritten: exp1_typeI.csv exp1_thresholds.csv exp1_F0.rds exp1_null2.csv exp1_qq.png\n")
