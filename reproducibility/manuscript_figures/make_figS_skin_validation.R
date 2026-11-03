## Supplementary figure: skin application, drawn from the values reported by Skin_MC903_MERFISH/code/skin_main_results.Rmd
## and skin_application_report.html (rendered reports; the Rmd sources recompute them from output/morphogam_auto). Descriptive only:
## one section per time point (atlas Mouse_ID = strip label; no within-time-point replication).
suppressMessages({library(dplyr); library(readr); library(ggplot2); library(patchwork)})
OUT <- "."   # run from this folder: Rscript make_figS_skin_validation.R
corr <- tribble(~reference, ~day, ~S, ~rank,
  "Papillary", 3, 0.224, 5,  "Papillary", 6, 0.463, 1,  "Papillary", 12, 0.507, 1,
  "Reticular", 3, 0.127, 6,  "Reticular", 6, 0.385, 1,  "Reticular", 12, 0.409, 1,
  "Subcutis", 6, 0.159, 5,   "Subcutis", 12, 0.179, 3,  "Subfascia", 6, 0.063, 6, "Subfascia", 12, 0.083, 3)
genes <- tribble(~gene, ~layer, ~strip, ~day, ~run, ~cells, ~frac_expressing, ~delta, ~se, ~z,
  "Cxcl12", "PAP", "d6_MC903_neck", 6, "20241007", 619, 0.32, 0.036, 0.012, 3.1,
  "Cxcl12", "PAP", "d9_MC903_Neck", 9, "20240408", 848, 0.25, 0.033, 0.011, 3.0,
  "Cxcl12", "PAP", "d12_MC903_neck", 12, "20241007", 208, 0.13, 0.073, 0.030, 2.4,
  "Cxcl12", "RET", "d6_MC903_neck", 6, "20241007", 377, 0.34, -0.016, 0.013, -1.2,
  "Cxcl12", "RET", "d9_MC903_Neck", 9, "20240408", 775, 0.34, 0.006, 0.006, 1.0,
  "Cxcl12", "RET", "d12_MC903_neck", 12, "20241007", 106, 0.26, 0.054, 0.040, 1.3,
  "Apcdd1", "PAP", "d3_MC903_neck", 3, "20241007", 49, 0.59, 0.009, 0.029, 0.3,
  "Apcdd1", "PAP", "d6_MC903_neck", 6, "20241007", 619, 0.16, -0.036, 0.013, -2.7,
  "Apcdd1", "PAP", "d12_MC903_neck", 12, "20241007", 208, 0.26, -0.042, 0.026, -1.6)
write_csv(corr, file.path(OUT, "skin_correspondence_values.csv")); write_csv(genes, file.path(OUT, "skin_gene_values.csv"))
PAL <- c(Papillary = "#e41a1c", Reticular = "#377eb8", Subcutis = "#e3b505", Subfascia = "#984ea3")
pA <- ggplot(corr, aes(day, S, colour = reference, group = reference)) + geom_line(linewidth = .6) + geom_point(size = 2) +
  geom_text(data = filter(corr, rank == 1), aes(label = "1"), vjust = -0.9, size = 2.3, show.legend = FALSE) +
  scale_colour_manual(values = PAL, name = "Healthy (day 0) reference") + scale_x_continuous(breaks = c(3, 6, 12)) + ylim(0, 0.6) +
  labs(x = "MC903 day of the FIB_ProInf section", y = "Positional similarity S(reference, FIB_ProInf)", title = "A  FIB_ProInf vs healthy fibroblast populations") + theme_bw(8) + theme(legend.position = "bottom", legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5))
gd <- genes %>% filter(gene == "Cxcl12") %>% mutate(layer = factor(layer, c("PAP", "RET"), c("papillary compartment", "reticular compartment")))
pB <- ggplot(gd, aes(day, delta, colour = layer)) + geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(ymin = delta - 1.96 * se, ymax = delta + 1.96 * se), position = position_dodge(width = 1.2), size = .35) +
  scale_colour_manual(values = c("#e41a1c", "#377eb8"), name = NULL) + scale_x_continuous(breaks = c(6, 9, 12)) +
  labs(x = "MC903 day (one section each)", y = expression(delta~"of"~italic(Cxcl12)~"within FIB_ProInf (±1.96 SE)"), title = "B  Position of Cxcl12 within FIB_ProInf") + theme_bw(8) + theme(legend.position = "bottom", legend.text = element_text(size = 6.5), legend.title = element_text(size = 6.5))
fig <- pA | pB
ggsave(file.path(OUT, "skin_validation_supplement.pdf"), fig, width = 180, height = 85, units = "mm", device = cairo_pdf)
ggsave(file.path(OUT, "skin_validation_supplement.png"), fig, width = 180, height = 85, units = "mm", dpi = 300, bg = "white"); cat("written\n")
