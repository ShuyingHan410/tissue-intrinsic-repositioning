## Figure 4 - simulation evaluation: A Type I error, B power, C specificity (standard null / population movement only /
## abundance only / repositioning only), D realized FDR and true-positive rate under BH.
## Every number comes from the frozen simulation result tables; nothing is simulated here.
##   Rscript make_fig4_simulation.R            (run from manuscript_figures/)
## Inputs: $SIM_DIR/exp1_typeI.csv, batch2_power.csv, batch2_spec.csv, batch2_fdr.csv (v16, summarise_*.R)
##         $SIM_DIR_V17/r8_absolute_centroid_baseline.csv (v17_revision/r8_absolute_centroid_baseline.R)
## Outputs next to this script: fig4_simulation.pdf/.png and fig4C_specificity_values.csv.
suppressMessages({library(dplyr); library(readr); library(ggplot2); library(patchwork); library(tidyr)})
V16 <- Sys.getenv("SIM_DIR", "../gene_spatial_cluster/simulation/v16"); V17 <- Sys.getenv("SIM_DIR_V17", "../gene_spatial_cluster/simulation/v17_revision"); OUT <- "."
TEAL <- "#1b7f79"; TEAL_L <- "#8fc9c4"; INK <- "#333333"; NAVY <- "#1f3b73"; TAUCOL <- c(`0.05` = "#9ed3cf", `0.15` = "#3a9f98", `0.30` = "#0f4f4b")
ci <- function(p, n = 200) 1.96 * sqrt(p * (1 - p) / n); pctlab <- scales::label_percent(accuracy = 1)
## House style shared with Figures 2, 3 and 5 (matplotlib scripts in this folder): Helvetica; bold 11 pt panel letter at the
## top-left corner; centred navy bold ~7 pt panel title on the same line; plain dark sub-panel labels; 5.5-6 pt axis text.
thm <- theme_minimal(base_size = 7, base_family = "Helvetica") + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(colour = "grey90", linewidth = .3), axis.line.x = element_line(colour = "grey55", linewidth = .35),
  axis.ticks = element_line(colour = "grey55", linewidth = .35), axis.ticks.length = unit(1.5, "pt"), axis.text = element_text(colour = INK, size = 5.6),
  axis.title = element_text(colour = INK, size = 6), strip.text = element_text(colour = INK, size = 6, hjust = .5, margin = margin(b = 2)),
  legend.title = element_text(size = 5.8, colour = INK), legend.text = element_text(size = 5.4, colour = INK), legend.key.size = unit(7, "pt"),
  plot.title = element_text(face = "bold", colour = NAVY, size = 7.2, hjust = .5, margin = margin(t = 1, b = 5)), plot.title.position = "plot",
  plot.tag = element_text(face = "bold", size = 11, colour = "black"), plot.tag.location = "plot", plot.tag.position = c(0, 1),
  plot.margin = margin(5, 6, 3, 3))
t1 <- read_csv(file.path(V16, "exp1_typeI.csv"), show_col_types = FALSE) %>% filter(method == "proposed", null == "null1")
pw <- read_csv(file.path(V16, "batch2_power.csv"), show_col_types = FALSE) %>% filter(method == "proposed", k == 10, Mratio == 1)
fd <- read_csv(file.path(V16, "batch2_fdr.csv"), show_col_types = FALSE) %>% filter(method == "proposed", mode == "native")
sp <- read_csv(file.path(V16, "batch2_spec.csv"), show_col_types = FALSE) %>% filter(method == "proposed")
r8 <- read_csv(file.path(V17, "r8_absolute_centroid_baseline.csv"), show_col_types = FALSE)
## ---- A
dA <- t1 %>% mutate(k = factor(k, c(10, 20))); NEFF <- min(dA$n_eff)
pA <- ggplot(dA, aes(tau, alpha, colour = k, group = k)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = .05 - ci(.05, NEFF), ymax = .05 + ci(.05, NEFF), fill = "grey88", alpha = .6) +
  geom_hline(yintercept = .05, linetype = "22", colour = "grey45", linewidth = .4) + geom_line(linewidth = .6) + geom_pointrange(aes(ymin = lo, ymax = hi), size = .25, linewidth = .4) +
  scale_colour_manual(values = c(`10` = TEAL, `20` = TEAL_L), name = "Sections per condition") + scale_x_continuous(breaks = c(0, .05, .15, .30)) +
  scale_y_continuous(labels = pctlab, limits = c(0, .10), breaks = seq(0, .10, .02), expand = expansion(0)) +
  labs(x = expression("Between-section heterogeneity " * tau[theta]), y = "Type I error (p < 0.05)", title = "Type I error, independent sections", tag = "A") +
  thm + theme(legend.position = c(.98, .98), legend.justification = c(1, 1), legend.background = element_rect(fill = "white", colour = NA))
## ---- B
dB <- pw %>% filter(tau %in% c(.05, .15, .30)) %>% transmute(scen, D, tau, power = reported) %>% mutate(tau = factor(sprintf("%.2f", tau), names(TAUCOL)), scen = factor(scen, c("gained", "reversed"), c("Gained", "Reversed")))
pB <- ggplot(dB, aes(D, power, colour = tau, linetype = tau, shape = tau, group = tau)) + geom_line(linewidth = .6) +
  geom_pointrange(aes(ymin = pmax(0, power - ci(power)), ymax = pmin(1, power + ci(power))), size = .25, linewidth = .4, linetype = "solid") + facet_wrap(~scen, nrow = 1) +
  scale_colour_manual(values = TAUCOL, name = expression(tau[theta])) + scale_linetype_manual(values = c(`0.05` = "solid", `0.15` = "42", `0.30` = "11"), name = expression(tau[theta])) +
  scale_shape_manual(values = c(`0.05` = 16, `0.15` = 17, `0.30` = 15), name = expression(tau[theta])) +
  scale_x_continuous(breaks = c(.025, .05, .075), labels = c("0.025", "0.050", "0.075"), limits = c(.02, .08)) + scale_y_continuous(labels = pctlab, limits = c(0, 1.02), breaks = seq(0, 1, .25), expand = expansion(0)) +
  labs(x = expression("Displacement |" * Delta * "|"), y = "Power (p < 0.05)", title = "Detection power, 10 sections per condition", tag = "B") + thm + theme(legend.position = "right")
## ---- C  specificity: standard null (Exp 1, tau .15, k 10), movement-only (R8 null3A pair), abundance-only (Exp 3 de_only), repositioning-only (Exp 3 repo_only)
null15 <- dA %>% filter(k == "10", tau == 0.15)
mv <- r8 %>% filter(scen == "null3A"); stopifnot(nrow(mv) == 1); mv$nominal <- mv$proposed_rej / 100
de <- sp %>% filter(scen == "de_only"); rp <- sp %>% filter(scen == "repo_only")
getp <- function(d) { cn <- intersect(c("nominal", "reported", "alpha", "rate", "rej"), names(d)); d[[cn[1]]] }
dC <- tibble(scenario = factor(c("Standard\nnull", "Population\nmoves +0.20", "Abundance\ndoubles", "Repositioning\n|Δ| = 0.05"),
                               levels = c("Standard\nnull", "Population\nmoves +0.20", "Abundance\ndoubles", "Repositioning\n|Δ| = 0.05")),
             rate = c(null15$alpha, getp(mv), getp(de), getp(rp)), n = c(NA, 200, 200, 200), kind = c("null", "null", "null", "alternative"))
dC <- dC %>% mutate(lo = ifelse(is.na(n), null15$lo, pmax(0, rate - ci(rate, n))), hi = ifelse(is.na(n), null15$hi, pmin(1, rate + ci(rate, n))))
write_csv(dC %>% mutate(source = c("v16 Exp 1 null1 tau 0.15 k 10 (50 x 200)", "v17 R8 null3A (200 reps)", "v16 Exp 3 de_only (200 reps)", "v16 Exp 3 repo_only (200 reps)")), file.path(OUT, "fig4C_specificity_values.csv"))
pC <- ggplot(dC, aes(scenario, rate, fill = kind)) + geom_hline(yintercept = .05, linetype = "22", colour = "grey45", linewidth = .4) + geom_col(width = .6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .18, linewidth = .4, colour = INK) + geom_text(aes(label = pctlab(rate), y = hi + .045), size = 1.9, colour = INK, family = "Helvetica") +
  scale_fill_manual(values = c(null = "grey70", alternative = TEAL), guide = "none") + scale_y_continuous(labels = pctlab, limits = c(0, 1.10), breaks = seq(0, 1, .25), expand = expansion(0)) +
  labs(x = "Simulated scenario (relative gene position unchanged in the nulls)", y = "Rejection rate (p < 0.05)", title = "Specificity to population-relative displacement", tag = "C") + thm + theme(axis.text.x = element_text(size = 5.6, lineheight = .9))
## ---- D  FDR / TPR
dD <- fd %>% transmute(k = factor(k), FDR = FDR_q05, TPR = TPR_q05, se = mc_se_q05) %>% pivot_longer(c(FDR, TPR), names_to = "stat", values_to = "v") %>%
  mutate(lo = ifelse(stat == "FDR", pmax(0, v - 1.96 * se), NA), hi = ifelse(stat == "FDR", v + 1.96 * se, NA),
         ylab = ifelse(stat == "FDR", hi + 0.004, v + 0.035),
         stat = factor(stat, c("FDR", "TPR"), c("Realized FDR", "True-positive rate")))
pD <- ggplot(dD, aes(k, v, fill = k)) + geom_col(width = .55) + geom_errorbar(aes(ymin = lo, ymax = hi), width = .15, linewidth = .4, colour = INK, na.rm = TRUE) +
  geom_hline(data = tibble(stat = factor("Realized FDR", levels(dD$stat)), y = .05), aes(yintercept = y), linetype = "22", colour = "grey45", linewidth = .4) +
  geom_text(aes(label = scales::label_percent(accuracy = 0.1)(v), y = ylab), size = 1.9, colour = INK, family = "Helvetica") + facet_wrap(~stat, scales = "free_y") +
  scale_fill_manual(values = c(`10` = TEAL, `20` = TEAL_L), guide = "none") + scale_y_continuous(labels = pctlab, expand = expansion(mult = c(0, .08))) +
  labs(x = "Sections per condition", y = "Rate", title = "Multiple testing, BH at q = 0.05", tag = "D") + thm
fig <- (pA | pB) / (pC | pD) + plot_layout(widths = c(1, 1.2), heights = c(1, 1))
## 165 x 112 mm: included at full text width (6.5 in) so every label keeps its design size (>= 6.5 pt).
ggsave(file.path(OUT, "fig4_simulation.pdf"), fig, width = 165, height = 112, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "fig4_simulation.png"), fig, width = 165, height = 112, units = "mm", dpi = 300, bg = "white")
print(dC); print(dD); cat("written\n")
