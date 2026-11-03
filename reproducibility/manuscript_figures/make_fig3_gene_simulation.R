## Figure 3 (main, 2 panels) + two supplementary figures for the gene repositioning simulations.
## Every value is read from the frozen v16 result files (simulation/v16/*.csv) and, for the absolute-centroid
## baseline only, from the v17 R8 run. Nothing is simulated, smoothed or interpolated here.
##   MAIN   fig3_gene_simulation_main.pdf/.png    A Type I error vs tau_theta (Exp 1, k 10/20)   B power vs |Delta| (Exp 2)
##   SUPP   figS_specificity.pdf/.png             three paired comparisons (R8 null3A; Exp 3 de_only; Exp 3 repo_only)
##   SUPP   figS_fdr_ablation.pdf/.png            A native-BH FDR (Exp 6)   B REML+KH vs fixed-effect pooling (Exp 4B)
## Usage: Rscript make_fig3_gene_simulation.R
suppressMessages({library(dplyr); library(readr); library(ggplot2); library(patchwork); library(tidyr)})
V16 <- "../gene_spatial_cluster/simulation/v16"; V17 <- "../gene_spatial_cluster/simulation/v17_revision"
TEAL <- "#1b7f79"; TEAL_L <- "#8fc9c4"; GREY <- "#8c8c8c"; INK <- "#333333"
TAUCOL <- c(`0.05` = "#9ed3cf", `0.15` = "#3a9f98", `0.30` = "#0f4f4b")
ci <- function(p, n = 200) 1.96 * sqrt(p * (1 - p) / n)
pctlab <- scales::label_percent(accuracy = 1)
thm <- theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = .3),
        axis.line.x = element_line(colour = "grey55", linewidth = .35), axis.ticks = element_line(colour = "grey55", linewidth = .35),
        axis.text = element_text(colour = INK, size = 7), axis.title = element_text(colour = INK, size = 7.5),
        strip.text = element_text(face = "bold", colour = INK, size = 8, hjust = 0),
        legend.title = element_text(size = 7), legend.text = element_text(size = 7), legend.key.size = unit(9, "pt"),
        plot.title = element_text(face = "bold", colour = INK, size = 9, hjust = 0),
        plot.tag = element_text(face = "bold", size = 12), plot.margin = margin(4, 8, 2, 2))

## ---------- frozen data ----------
t1  <- read_csv(file.path(V16, "exp1_typeI.csv"), show_col_types = FALSE) %>% filter(method == "proposed", null == "null1")
abB <- read_csv(file.path(V16, "ablB_k10.csv"), show_col_types = FALSE)
pw  <- read_csv(file.path(V16, "batch2_power.csv"), show_col_types = FALSE) %>% filter(method == "proposed", k == 10)
fd  <- read_csv(file.path(V16, "batch2_fdr.csv"), show_col_types = FALSE) %>% filter(method == "proposed", mode == "native")
sp  <- read_csv(file.path(V16, "batch2_spec.csv"), show_col_types = FALSE) %>% filter(method %in% c("proposed", "pbulk_de"))
r8  <- read_csv(file.path(V17, "r8_absolute_centroid_baseline.csv"), show_col_types = FALSE)

## ================= MAIN FIGURE =================
dA <- t1 %>% mutate(k = factor(k, c(10, 20), c("10", "20")))
## reference band: Monte Carlo 95% interval of the rejection proportion under a true 5% rate, using the smallest
## cluster-robust effective number of tests among the plotted cells (dependence of the 50 genes within a replicate)
NEFF <- min(dA$n_eff)
pA <- ggplot(dA, aes(tau, alpha, colour = k, group = k)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = .05 - ci(.05, NEFF), ymax = .05 + ci(.05, NEFF), fill = "grey88", alpha = .6) +
  geom_hline(yintercept = .05, linetype = "22", colour = "grey45", linewidth = .4) +
  annotate("text", x = .30, y = .05 + ci(.05, NEFF), label = "95% Monte Carlo reference\ninterval under the nominal null", hjust = 1, vjust = -.3, size = 2.1, lineheight = .85, colour = "grey40") +
  geom_line(linewidth = .7) + geom_pointrange(aes(ymin = lo, ymax = hi), size = .3, linewidth = .4) +
  scale_colour_manual(values = c(`10` = TEAL, `20` = TEAL_L), name = "Sections per condition") +
  scale_x_continuous(breaks = c(0, .05, .15, .30)) +
  scale_y_continuous(labels = pctlab, limits = c(0, .10), breaks = seq(0, .10, .02), expand = expansion(0)) +
  labs(x = expression("Between-section heterogeneity, " * tau[theta]), y = "Empirical Type I error (p < 0.05)", title = "Empirical Type I error", tag = "A") +
  thm + theme(legend.position = c(.98, .98), legend.justification = c(1, 1), legend.direction = "vertical",
              legend.background = element_rect(fill = "white", colour = NA))

dB <- pw %>% filter(tau %in% c(.05, .15, .30)) %>% transmute(scen, D, tau, power = reported) %>%
  mutate(tau = factor(sprintf("%.2f", tau), names(TAUCOL)),
         scen = factor(scen, c("gained", "reversed"), c("Gained localization", "Reversed localization")))
pB <- ggplot(dB, aes(D, power, colour = tau, linetype = tau, shape = tau, group = tau)) +
  geom_line(linewidth = .7) +
  geom_pointrange(aes(ymin = pmax(0, power - ci(power)), ymax = pmin(1, power + ci(power))), size = .3, linewidth = .4, linetype = "solid") +
  facet_wrap(~scen, nrow = 1) +
  scale_colour_manual(values = TAUCOL, name = expression(tau[theta])) +
  scale_linetype_manual(values = c(`0.05` = "solid", `0.15` = "42", `0.30` = "11"), name = expression(tau[theta])) +
  scale_shape_manual(values = c(`0.05` = 16, `0.15` = 17, `0.30` = 15), name = expression(tau[theta])) +
  scale_x_continuous(breaks = c(.025, .05, .075), labels = c("0.025", "0.050", "0.075"), limits = c(.02, .08)) +
  scale_y_continuous(labels = pctlab, limits = c(0, 1.02), breaks = seq(0, 1, .25), expand = expansion(0)) +
  labs(x = expression("Displacement magnitude, |" * Delta * "|"), y = "Empirical power (nominal p < 0.05)", title = "Detection power", tag = "B") +
  thm + theme(legend.position = "right", legend.justification = "top")
main <- pA + pB + plot_layout(widths = c(.35, .65))
ggsave("legacy_fig3_gene_simulation_main.pdf", main, width = 180, height = 72, units = "mm", device = cairo_pdf)
ggsave("legacy_fig3_gene_simulation_main.png", main, width = 180, height = 72, units = "mm", dpi = 300, bg = "white")

## ================= SUPP: specificity =================
get <- function(sc, m) sp %>% filter(scen == sc, method == m) %>% pull(nominal)
dD <- bind_rows(
  tibble(panel = "population moves (+0.20),\nrelative gene position unchanged", stat = c("ours (δ)", "absolute transcript\ncentroid"),
         rej = c(r8$proposed_rej[r8$scen == "null3A"], r8$absolute_rej[r8$scen == "null3A"]) / 100),
  tibble(panel = "gene abundance doubles,\nno movement", stat = c("ours (δ)", "pseudobulk DE"), rej = c(get("de_only", "proposed"), get("de_only", "pbulk_de"))),
  tibble(panel = "directional repositioning only,\n|Δ| = 0.05", stat = c("ours (δ)", "pseudobulk DE"), rej = c(get("repo_only", "proposed"), get("repo_only", "pbulk_de")))) %>%
  mutate(panel = factor(panel, unique(panel)), ours = stat == "ours (δ)", stat = factor(stat, c("ours (δ)", "absolute transcript\ncentroid", "pseudobulk DE")))
pD <- ggplot(dD, aes(stat, rej, fill = ours)) + geom_col(width = .6) +
  geom_errorbar(aes(ymin = pmax(0, rej - ci(rej)), ymax = pmin(1, rej + ci(rej))), width = .15, linewidth = .35, colour = INK) +
  geom_hline(yintercept = .05, linetype = "22", colour = "grey45", linewidth = .4) +
  geom_text(aes(y = pmin(1, rej + ci(rej)), label = sprintf("%.1f%%", 100 * rej)), vjust = -.5, size = 2.4, colour = INK) +
  facet_wrap(~panel, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = TEAL, `FALSE` = GREY), guide = "none") +
  scale_y_continuous(labels = pctlab, limits = c(0, 1.1), breaks = seq(0, 1, .25), expand = expansion(0)) +
  labs(x = NULL, y = "Rejection rate (p < 0.05)") + thm + theme(strip.text = element_text(size = 7.5, lineheight = .9), axis.text.x = element_text(lineheight = .85))
ggsave("legacy_figS_specificity.pdf", pD, width = 180, height = 70, units = "mm", device = cairo_pdf)
ggsave("legacy_figS_specificity.png", pD, width = 180, height = 70, units = "mm", dpi = 300, bg = "white")

## ================= SUPP: FDR + ablation =================
dC <- fd %>% transmute(k = factor(k, c(10, 20), c("10", "20")), FDR = FDR_q05, se = mc_se_q05, TPR = TPR_q05)
pC <- ggplot(dC, aes(k, FDR)) + geom_col(width = .55, fill = TEAL) +
  geom_errorbar(aes(ymin = FDR - 1.96 * se, ymax = FDR + 1.96 * se), width = .15, linewidth = .35, colour = INK) +
  geom_hline(yintercept = .05, linetype = "22", colour = "grey45", linewidth = .4) +
  geom_text(aes(y = FDR + 1.96 * se, label = sprintf("%.1f%%", 100 * FDR)), vjust = -.5, size = 2.4, colour = INK) +
  geom_text(aes(y = .003, label = sprintf("TPR\n%.0f%%", 100 * TPR)), vjust = 0, size = 2.1, lineheight = .85, colour = "white") +
  scale_y_continuous(labels = pctlab, limits = c(0, .07), breaks = seq(0, .06, .02), expand = expansion(0)) +
  labs(x = "Sections per condition", y = "Realized FDR (BH, q = 0.05)", title = "FDR control", tag = "A") + thm
dE <- abB %>% mutate(pooling = factor(pooling, c("REML+KH", "fixed"), c("REML + Knapp–Hartung", "fixed-effect")),
                     tau = factor(tau, c(0, .18), c("τθ = 0", "τθ = 0.18 (I² ≈ 60%)")), what = ifelse(D == 0, "null", sprintf("|Δ| = %.3f", D)),
                     what = factor(what, c("null", "|Δ| = 0.025", "|Δ| = 0.050", "|Δ| = 0.075")))
pE <- ggplot(dE, aes(what, reject, fill = pooling)) + geom_col(position = position_dodge(.7), width = .62) +
  geom_errorbar(aes(ymin = pmax(0, reject - ci(reject)), ymax = pmin(1, reject + ci(reject))), position = position_dodge(.7), width = .15, linewidth = .35, colour = INK) +
  geom_hline(yintercept = .05, linetype = "22", colour = "grey45", linewidth = .4) +
  facet_wrap(~tau, nrow = 1) + scale_fill_manual(values = c(TEAL, GREY), name = NULL) +
  scale_y_continuous(labels = pctlab, limits = c(0, 1.05), breaks = seq(0, 1, .25), expand = expansion(0)) +
  labs(x = NULL, y = "Rejection rate (p < 0.05)", title = "Pooling ablation", tag = "B") +
  thm + theme(legend.position = "bottom", legend.margin = margin(-4, 0, 0, 0))
supp2 <- pC + pE + plot_layout(widths = c(.5, 1.5))
ggsave("figS_fdr_ablation.pdf", supp2, width = 180, height = 72, units = "mm", device = cairo_pdf)
ggsave("figS_fdr_ablation.png", supp2, width = 180, height = 72, units = "mm", dpi = 300, bg = "white")

## ---------- numbers behind the figures ----------
out <- bind_rows(
  dA %>% transmute(figure = "Fig3A", x = tau, group = paste("k", k), value = alpha, lo, hi, n_rep = R),
  dB %>% transmute(figure = "Fig3B", x = D, group = paste(scen, "tau", tau), value = power, lo = pmax(0, power - ci(power)), hi = pmin(1, power + ci(power)), n_rep = 200),
  dD %>% transmute(figure = "FigS_specificity", x = NA, group = paste(gsub("\n", " ", panel), gsub("\n", " ", stat), sep = " :: "), value = rej, lo = pmax(0, rej - ci(rej)), hi = pmin(1, rej + ci(rej)), n_rep = 200),
  dC %>% transmute(figure = "FigS_fdr", x = NA, group = paste("k", k), value = FDR, lo = FDR - 1.96 * se, hi = FDR + 1.96 * se, n_rep = 50),
  dE %>% transmute(figure = "FigS_ablation", x = D, group = paste(tau, pooling, what, sep = " :: "), value = reject, lo = pmax(0, reject - ci(reject)), hi = pmin(1, reject + ci(reject)), n_rep = 200))
dir.create("fig23_assets", showWarnings = FALSE); write_csv(out, "fig23_assets/fig3_gene_sim_data.csv")
cat("written figS_fdr_ablation (submitted supplement) + legacy_fig3_gene_simulation_main, legacy_figS_specificity (not used) (.pdf/.png) + fig23_assets/fig3_gene_sim_data.csv\n")
