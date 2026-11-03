## Summarise batch 2 (Experiments 2, 3, 6) under the v16 §0.5 rules, using the design-specific
## null distributions from Experiment 1 / 1b (exp1_thresholds.csv, exp1_F0.rds).
## Produces:
##   batch2_power.csv    Exp 2: nominal and matched power per method x cell, testability columns
##   batch2_spec.csv     Exp 3: rejection rates per method x scenario (labelled, never "Type I")
##   batch2_fdr.csv      Exp 6: native BH and F0-calibrated BH, realised FDR / TPR at q in {.05,.10}; oracle q*
##   batch2_fig_*.png    power vs |D| by tau; nominal vs matched at |D| = .05; FDR curves
## Method roles (§0.3): primary = proposed, gam_fs, despace2; pbulk_de = control; gam_naive = Type I panels only.
suppressMessages({source("engine_v16.R"); library(ggplot2)})
stopifnot(file.exists("exp1_thresholds.csv"), file.exists("exp1_F0.rds"))
thr <- read_csv("exp1_thresholds.csv", show_col_types = FALSE)
F0  <- readRDS("exp1_F0.rds")
t1  <- read_csv("exp1_typeI.csv", show_col_types = FALSE) %>% filter(null == "null1")
cut_of <- function(m, k, tau) { z <- thr %>% filter(method == m, k == !!k, abs(tau - !!tau) < 1e-9); if (nrow(z)) z$c_0.05 else NA_real_ }
nominal_ok <- function(m, k, tau) { z <- t1 %>% filter(method == m, k == !!k, abs(tau - !!tau) < 1e-9); if (nrow(z)) isTRUE(z$pass_3_7) else NA }
calib_p <- function(p, m, k, tau) { key <- sprintf("%s|k%d|tau%.2f", m, k, tau); f <- F0[[key]]; if (is.null(f)) return(rep(NA_real_, length(p)))
  ecdf(f)(p) }                                                    # p~ = F0(p), uniform under a valid null
ci <- function(a, n) 1.96 * sqrt(a * (1 - a) / n)
PRIMARY <- c("proposed", "gam_fs", "despace2"); ALL5 <- c(PRIMARY, "pbulk_de", "gam_naive")

cells <- lapply(list.files("batch2_out", "^cell_\\d+[.]rds$", full.names = TRUE), readRDS)
## ---- Experiments 2 and 3 --------------------------------------------------
rows <- list()
for (z in cells) { if (z$cell$exp == "6") next
  P <- z$P; n <- nrow(P)
  for (m in ALL5) {
    p <- P[, m]; cc <- cut_of(m, z$cell$k, z$cell$tau)
    rows[[length(rows) + 1]] <- data.frame(z$cell, method = m, n = n,
      nominal = mean(p < 0.05, na.rm = TRUE), matched = if (is.na(cc)) NA else mean(p < cc, na.rm = TRUE), cutoff = cc,
      nominal_ok = nominal_ok(m, z$cell$k, z$cell$tau),
      P_testable = if (m == "proposed") mean(P[, "k0"] >= 3 & P[, "k9"] >= 3 & is.finite(p)) else mean(is.finite(p)),
      k0_med = if (m == "proposed") median(P[, "k0"]) else NA, k9_med = if (m == "proposed") median(P[, "k9"]) else NA,
      I2_realised = if ("I2" %in% colnames(P)) median(P[, "I2"], na.rm = TRUE) else NA,
      aux = if ("aux" %in% colnames(P)) median(P[, "aux"], na.rm = TRUE) else NA)
  }
}
tab <- bind_rows(rows) %>% mutate(reported = ifelse(nominal_ok %in% TRUE, nominal, matched),
                                  rule = ifelse(nominal_ok %in% TRUE, "nominal p<0.05", sprintf("matched p<%.4f", cutoff)),
                                  P_correct_discovery = P_testable * reported)
pw <- tab %>% filter(exp == "2"); sp <- tab %>% filter(exp == "3")
write_csv(pw, "batch2_power.csv"); write_csv(sp, "batch2_spec.csv")
options(width = 220)
cat("\n## Experiment 2 — power (reported rule per §0.5), |D| x tau, primary methods\n")
print(as.data.frame(pw %>% filter(method %in% PRIMARY) %>% mutate(v = sprintf("%.3f %s", reported, ifelse(nominal_ok %in% TRUE, "", "(m)"))) %>%
  select(scen, D, tau, method, v) %>% pivot_wider(names_from = method, values_from = v) %>% arrange(scen, D, tau)))
cat("   (m) = reported at the design-specific matched threshold because the method was non-nominal there in Experiment 1\n")
cat("\n## Experiment 3 — rejection rates (NOT Type I for GAM/DESpace2 on nd; NOT Type I for pbulk on de_only)\n")
print(as.data.frame(sp %>% select(scen, method, nominal, matched, P_testable, aux) %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

## ---- Experiment 6 — FDR ---------------------------------------------------
fdr_rows <- list(); curves <- list()
for (z in cells) { if (z$cell$exp != "6") next
  k <- z$cell$k; tau <- z$cell$tau; truth <- z$truth
  for (m in c(PRIMARY, "pbulk_de")) for (mode in c("native", "calibrated")) {
    A <- sapply(z$P, function(P) { p <- P[, m]; if (mode == "calibrated") p <- calib_p(p, m, k, tau); p.adjust(p, "BH") })
    if (all(is.na(A))) next
    ev <- function(q) { fdp <- apply(A, 2, function(a) { s <- a < q; R <- sum(s, na.rm = TRUE); if (R > 0) sum(s & !truth, na.rm = TRUE) / R else 0 })
      tpr <- apply(A, 2, function(a) sum(a < q & truth, na.rm = TRUE) / sum(truth)); c(FDR = mean(fdp), TPR = mean(tpr), se = sd(fdp) / sqrt(length(fdp))) }
    a5 <- ev(.05); a10 <- ev(.10)
    fdr_rows[[length(fdr_rows) + 1]] <- data.frame(k = k, tau = tau, method = m, mode = mode, n_rep = length(z$P),
      FDR_q05 = a5["FDR"], TPR_q05 = a5["TPR"], FDR_q10 = a10["FDR"], TPR_q10 = a10["TPR"], mc_se_q05 = a5["se"])
    qs <- seq(.002, .2, by = .002); cv <- t(vapply(qs, function(q) ev(q)[1:2], numeric(2)))
    curves[[length(curves) + 1]] <- data.frame(k = k, method = m, mode = mode, q = qs, FDR = cv[, 1], TPR = cv[, 2])
  }
}
if (length(fdr_rows)) {
  fd <- bind_rows(fdr_rows); cv <- bind_rows(curves)
  orc <- cv %>% filter(mode == "native", method %in% PRIMARY) %>% group_by(k, method) %>% slice_min(abs(FDR - .05), n = 1) %>%
    transmute(k, method, q_star = q, FDR_at_qstar = FDR, TPR_at_qstar = TPR)
  fd <- fd %>% left_join(orc, by = c("k", "method")); write_csv(fd, "batch2_fdr.csv")
  cat("\n## Experiment 6 — realised FDR / TPR under BH: native vs F0-calibrated p (§0.5 step 3)\n")
  print(as.data.frame(fd %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
  g <- ggplot(cv %>% filter(method %in% PRIMARY), aes(q, FDR, colour = method, linetype = mode)) + geom_abline(linetype = "dotted", colour = "grey50") +
    geom_line() + facet_wrap(~k, labeller = label_both) + labs(x = "nominal BH q", y = "realised FDR") + theme_minimal(base_size = 10)
  ggsave("batch2_fig_fdr.png", g, width = 9, height = 4, dpi = 150)
}
## ---- figures for Exp 2 --------------------------------------------------------
g1 <- ggplot(pw %>% filter(method %in% PRIMARY), aes(D, reported, colour = method)) + geom_line() + geom_point() +
  facet_grid(scen ~ tau, labeller = label_both) + ylim(0, 1) + labs(x = "|D|", y = "power (reported rule)") + theme_minimal(base_size = 10)
ggsave("batch2_fig_power.png", g1, width = 11, height = 5.5, dpi = 150)
g2 <- pw %>% filter(method %in% PRIMARY, abs(D - .05) < 1e-9) %>% group_by(method, tau) %>% summarise(nominal = mean(nominal), matched = mean(matched, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(c(nominal, matched), names_to = "type", values_to = "power") %>%
  ggplot(aes(method, power, fill = method, alpha = type)) + geom_col(position = position_dodge(.7), width = .6) + facet_wrap(~tau, labeller = label_both, nrow = 1) +
  scale_alpha_manual(values = c(nominal = .45, matched = 1)) + ylim(0, 1) + labs(y = "power at |D| = 0.05 (gained/reversed mean)") + theme_minimal(base_size = 10) + theme(axis.text.x = element_blank())
ggsave("batch2_fig_matched.png", g2, width = 10, height = 3.8, dpi = 150)
cat("\nwritten: batch2_power.csv batch2_spec.csv batch2_fdr.csv batch2_fig_power.png batch2_fig_matched.png batch2_fig_fdr.png\n")
