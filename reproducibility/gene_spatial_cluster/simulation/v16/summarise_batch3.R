## Summarise batch 3 (Experiments 5, 7, 7b, 8, 9, 10, 11, 6b) under the v16 §0.5 rules.
## Thresholds: for pseudo-arm cells at (k, tau) present in Experiment 1/1b, the design-specific c_{m,k,tau};
## for other substrates / abundances the cell's OWN null1 partner in this batch (same exp, sub, tau, M_B, k)
## supplies the cutoff, per the "null threshold and power scenario must share the design" rule.
## Outputs: batch3_table.csv (every cell x method: nominal, matched, rule, testability, realised I2),
##          batch3_exp9.csv (lambda diagnostics), batch3_fdr6b.csv, batch3_fig_*.png
suppressMessages({source("engine_v16.R"); library(ggplot2)})
thr <- read_csv("exp1_thresholds.csv", show_col_types = FALSE); t1 <- read_csv("exp1_typeI.csv", show_col_types = FALSE) %>% filter(null == "null1")
F0 <- readRDS("exp1_F0.rds")
cells <- lapply(list.files("batch3_out", "^cell_\\d+[.]rds$", full.names = TRUE), readRDS)
info <- bind_rows(lapply(cells, function(z) data.frame(z$cell, key = z$key, n_rep = z$n_rep)))
P3 <- c("proposed", "gam_fs", "despace2"); ALL <- c(P3, "pbulk_de", "gam_naive", "proposed_lambda")

## null partner within batch 3: same exp/sub/tau/M_B/k/sigma/mkind with scen null1 (or null3B for tilt)
partner_null <- function(cl) {
  cand <- Filter(function(z) { c2 <- z$cell; c2$exp == cl$exp && c2$sub == cl$sub && abs(c2$tau - cl$tau) < 1e-9 && c2$M_B == cl$M_B &&
    c2$k == cl$k && c2$sigma == cl$sigma && c2$mkind == cl$mkind && c2$scen %in% c("null1", "null3B") }, cells)
  if (length(cand)) cand[[1]] else NULL
}
cut_for <- function(m, cl) {
  ## Exp 8 / Exp 11 / Exp 7 at M=64 would use Exp 1 tables; everything else uses its in-batch null partner
  if (cl$sub == "pseudo" && cl$M_B == 64 && cl$sigma == 0 && cl$shape == "exp") { z <- thr %>% filter(method == m, k == cl$k, abs(tau - cl$tau) < 1e-9); if (nrow(z)) return(c(cut = z$c_0.05, src = "exp1")) }
  pn <- partner_null(cl); if (is.null(pn) || !(m %in% colnames(pn$P))) return(c(cut = NA, src = "none"))
  p <- pn$P[, m]; c(cut = unname(quantile(p[is.finite(p)], .05)), src = "batch3_null")
}
rows <- list()
for (z in cells) { if (z$cell$exp == "6b") next
  P <- z$P; n <- nrow(P)
  for (m in intersect(ALL, colnames(P))) { p <- P[, m]; cc <- cut_for(if (m == "proposed_lambda") "proposed" else m, z$cell)
    rows[[length(rows) + 1]] <- data.frame(z$cell, method = m, n = n,
      nominal = mean(p < .05, na.rm = TRUE), matched = if (is.na(cc["cut"])) NA else mean(p < as.numeric(cc["cut"]), na.rm = TRUE),
      cutoff = as.numeric(cc["cut"]), cut_src = cc["src"],
      P_testable = mean(P[, "k0"] >= 3 & P[, "k9"] >= 3 & is.finite(p)), k0_med = median(P[, "k0"]), k9_med = median(P[, "k9"]),
      I2_realised = median(P[, "I2"], na.rm = TRUE)) } }
tab <- bind_rows(rows) %>% mutate(P_correct = P_testable * ifelse(scen %in% c("null1", "null3A", "null3B"), NA, nominal))
write_csv(tab, "batch3_table.csv")
options(width = 230)
show <- function(e, extra = NULL) { d <- tab %>% filter(exp == e, method %in% c(P3, "proposed_lambda", "pbulk_de")) %>%
  mutate(v = sprintf("%.3f/%.3f", nominal, matched)) %>% select(sub, scen, D, tau, M_B, k, shape, sigma, mkind, method, v, P_testable, I2_realised) %>%
  pivot_wider(names_from = method, values_from = v); cat(sprintf("\n## Experiment %s  (nominal / matched)\n", e)); print(as.data.frame(d)) }
for (e in c("5", "7", "7b", "8", "10", "11")) show(e)
## Exp 9: lambda diagnostics
e9 <- bind_rows(lapply(Filter(function(z) z$cell$exp == "9", cells), function(z) { P <- z$P
  data.frame(z$cell[, c("scen", "tau", "sigma", "mkind")], n = nrow(P), proposed = mean(P[, "proposed"] < .05), proposed_lambda = mean(P[, "proposed_lambda"] < .05),
    gam_fs = mean(P[, "gam_fs"] < .05, na.rm = TRUE), despace2 = mean(P[, "despace2"] < .05, na.rm = TRUE),
    lam_raw_med = median(P[, "lam_raw"], na.rm = TRUE), floor_fires = mean(P[, "lam_applied"] <= 1 + 1e-9), lam_applied_med = median(P[, "lam_applied"]),
    blank_fp = mean(P[, "blank_fp"]), blank_z_sd = median(P[, "blank_z_sd"], na.rm = TRUE), blank_q001 = median(P[, "blank_p_q001"], na.rm = TRUE)) }))
write_csv(e9, "batch3_exp9.csv"); cat("\n## Experiment 9 — deployed pipeline (lambda on) x mouse-effect distribution\n"); print(as.data.frame(e9 %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
## Exp 6b
z6 <- Filter(function(z) z$cell$exp == "6b", cells)
if (length(z6)) { z <- z6[[1]]; truth <- z$truth; out <- list()
  for (m in c(P3, "pbulk_de")) for (mode in c("native", "calibrated")) {
    A <- sapply(z$P, function(P) { p <- P[, m]; if (mode == "calibrated") { f <- F0[[sprintf("%s|k10|tau0.15", m)]]; if (!is.null(f)) p <- ecdf(f)(p) }; p.adjust(p, "BH") })
    ev <- function(q) { fdp <- apply(A, 2, function(a) { s <- a < q; R <- sum(s, na.rm = TRUE); if (R > 0) sum(s & !truth, na.rm = TRUE) / R else 0 }); tpr <- apply(A, 2, function(a) sum(a < q & truth, na.rm = TRUE) / sum(truth)); c(mean(fdp), mean(tpr)) }
    a5 <- ev(.05); a10 <- ev(.10); out[[length(out) + 1]] <- data.frame(method = m, mode = mode, FDR_q05 = a5[1], TPR_q05 = a5[2], FDR_q10 = a10[1], TPR_q10 = a10[2]) }
  f6 <- bind_rows(out); write_csv(f6, "batch3_fdr6b.csv"); cat("\n## Experiment 6b — FDR with block-correlated genes (rho = .7, blocks of 20); compare batch2_fdr.csv k = 10\n"); print(as.data.frame(f6 %>% mutate(across(where(is.numeric), ~round(.x, 3))))) }
## figures: Exp 7 power vs M at tau .05; Exp 9 Type I by sigma x kind
g7 <- tab %>% filter(exp == "7", scen == "gained", abs(D - .05) < 1e-9, method %in% P3) %>% bind_rows(read_csv("batch2_power.csv", show_col_types = FALSE) %>% filter(scen == "gained", abs(D - .05) < 1e-9, method %in% P3) %>% mutate(M_B = 64) %>% select(tau, M_B, method, nominal, matched)) %>%
  ggplot(aes(factor(M_B), matched, colour = method, group = method)) + geom_line() + geom_point() + facet_wrap(~tau, labeller = label_both, nrow = 1) + ylim(0, 1) + labs(x = "M (molecules per slice)", y = "matched power, gained |D| = .05") + theme_minimal(base_size = 10)
ggsave("batch3_fig_exp7.png", g7, width = 10, height = 3.6, dpi = 150)
g9 <- e9 %>% filter(scen == "null1") %>% pivot_longer(c(proposed, proposed_lambda, gam_fs, despace2), names_to = "method", values_to = "typeI") %>%
  ggplot(aes(mkind, typeI, colour = method, group = method)) + geom_hline(yintercept = .05, linetype = "dotted") + geom_point(size = 2) + geom_line() + facet_grid(tau ~ sigma, labeller = label_both) + labs(x = "mouse-effect distribution", y = "Type I (Null 1)") + theme_minimal(base_size = 10)
ggsave("batch3_fig_exp9.png", g9, width = 9, height = 5, dpi = 150)
cat("\nwritten: batch3_table.csv batch3_exp9.csv batch3_fdr6b.csv batch3_fig_exp7.png batch3_fig_exp9.png\n")
