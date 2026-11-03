## Step 4 of the stratified analysis: the gene-call table behind Figure 5.
##   Rscript 04_figure5_calls.R [contrastsC.csv] [out.csv]
## Takes the D1_full rows of the stratified contrast table (02d_stratified_merge.R), applies
## Benjamini-Hochberg to the raw Welch p-value of Delta = delta(Day 9) - delta(Day 0) within each
## correspondence x anatomical-layer family (eight families in the colitis analysis), and writes
## gene_calls_D1full_noLambda.csv, the input of manuscript_figures/make_fig5_repositioning_v2.py.
## No blank-barcode lambda is applied here (the lambda diagnostics are in lambdaC.csv); the
## manuscript reports these q-values. The endpoint classes (gained / lost / reversed / amplified /
## dampened) use the unwidened REML + Knapp-Hartung endpoint intervals (the *_raw columns written by
## 02d_stratified_merge.R), as in the manuscript. 892 rows; 113 at q < 0.05.
suppressMessages({library(dplyr); library(readr)})
a <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(a) >= 1) a[1] else "contrastsC.csv"
OUT <- if (length(a) >= 2) a[2] else "gene_calls_D1full_noLambda.csv"
con <- read_csv(IN, show_col_types = FALSE, progress = FALSE)
if ("ci_lo_parent_D0_raw" %in% names(con)) {           # use the unwidened endpoint intervals and p-values
  con <- con %>% mutate(ci_lo_parent_D0 = ci_lo_parent_D0_raw, ci_hi_parent_D0 = ci_hi_parent_D0_raw,
                        ci_lo_daughter_D9 = ci_lo_daughter_D9_raw, ci_hi_daughter_D9 = ci_hi_daughter_D9_raw,
                        p_parent_D0 = p_parent_D0_raw, p_daughter_D9 = p_daughter_D9_raw)
}
pos <- function(lo, hi) is.finite(lo) & is.finite(hi) & (lo > 0 | hi < 0)
con <- con %>% mutate(d0_positioned = pos(ci_lo_parent_D0, ci_hi_parent_D0), d9_positioned = pos(ci_lo_daughter_D9, ci_hi_daughter_D9),
  pattern = case_when(route != "D1_full" ~ NA_character_,
    !d0_positioned &  d9_positioned ~ "gained",
     d0_positioned & !d9_positioned ~ "lost",
     d0_positioned &  d9_positioned & sign(delta_parent_D0) != sign(delta_daughter_D9) ~ "reversed",
     d0_positioned &  d9_positioned & abs(delta_daughter_D9) > abs(delta_parent_D0) ~ "amplified",
     d0_positioned &  d9_positioned ~ "dampened", TRUE ~ "weak-both-ends"))
calls <- con %>% filter(route == "D1_full", !grepl("^Blank", gene)) %>%
  group_by(lineage, layer) %>% mutate(q_nolambda = p.adjust(d1_p_raw, "BH")) %>% ungroup() %>%
  transmute(gene, lineage, layer, delta_parent_D0, delta_daughter_D9, se_parent_D0, se_daughter_D9,
            ci_lo_parent_D0, ci_hi_parent_D0, ci_lo_daughter_D9, ci_hi_daughter_D9,
            d1, d1_se, d1_p_raw, p_raw = d1_p_raw, q_nolambda, pattern, direction,
            k_slices_parent_D0, k_slices_daughter_D9, cells_parent_D0, cells_daughter_D9) %>%
  arrange(lineage, layer, gene)
write_csv(calls, OUT)
cat(sprintf("%s: %d comparisons in %d families; %d at q < 0.05 (%d distinct genes)\n", OUT, nrow(calls),
            n_distinct(paste(calls$lineage, calls$layer)), sum(calls$q_nolambda < 0.05), n_distinct(calls$gene[calls$q_nolambda < 0.05])))
