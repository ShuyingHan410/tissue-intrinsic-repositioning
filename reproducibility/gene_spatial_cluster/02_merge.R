# Step 2 (one job): bind all chunk outputs, calibrate on the 48 Blank probes, classify.
# Route-aware: D1_full / daughter_only / parent_only are judged against DIFFERENT nulls.
#   -> gene_spatial_lineage_results.csv   every gene x lineage x design (nothing dropped)
#   -> gene_spatial_lineage_pairs.csv     every slice pair (raw)
#   -> gene_spatial_lineage_slices.csv    every slice delta (raw)
#   -> hits.csv                           calibrated + classified candidates
OUT <- Sys.getenv("OUT", ".")
suppressMessages({library(dplyr); library(readr); library(tidyr)})
RES <- file.path(OUT, "results")
bind_parts <- function(pat) {
  fs <- list.files(RES, pat, full.names = TRUE)
  if (!length(fs)) return(tibble())
  bind_rows(lapply(fs, read_csv, show_col_types = FALSE, progress = FALSE))
}
res <- bind_parts("^part_.*[.]csv$")
prs <- bind_parts("^pairs_.*[.]csv$")
slc <- bind_parts("^slices_.*[.]csv$")
gi  <- read_csv(file.path(OUT, "genes_index.csv"), show_col_types = FALSE)
res <- res %>% left_join(gi %>% select(gene, is_blank), by = "gene")
write_csv(res, file.path(OUT, "gene_spatial_lineage_results.csv"))
write_csv(prs, file.path(OUT, "gene_spatial_lineage_pairs.csv"))
write_csv(slc, file.path(OUT, "gene_spatial_lineage_slices.csv"))
cat(sprintf("bound: %d summary | %d pair | %d slice rows\n", nrow(res), nrow(prs), nrow(slc)))
cat("\nroutes:\n"); print(table(res$route, res$lineage, useNA = "ifany"))

# ---- NULL 1: shift-based, for D1_full rows -------------------------------------------
# Blank barcodes carry no transcript, so their shift/consistency is pure noise at this
# tissue, depth and cell count. Thresholds are the 95th percentile OF THE BLANKS, computed
# PER LINEAGE (cell counts and depth differ between lineages).
thr_d1 <- res %>% filter(is_blank, design == "D1_lineage", is.finite(median_shift)) %>%
  group_by(lineage) %>%
  summarise(n_blank = n(),
            thr_cons  = quantile(consistency, .95, na.rm = TRUE),
            thr_shift = quantile(abs(median_shift), .95, na.rm = TRUE), .groups = "drop")
cat("\nblank null for D1 (per lineage):\n"); print(as.data.frame(thr_d1), row.names = FALSE)
if (!nrow(thr_d1)) stop(
  "BLANK NULL IS EMPTY - no rows with is_blank=TRUE and design='D1_lineage'.\n",
  "  Without it every gene fails clears_null and the run silently returns 0 hits.\n",
  "  Check: (a) genes_index.csv has is_blank=TRUE for the 48 Blank-NN columns,\n",
  "         (b) the chunks actually contain the Blank-NN columns (they are analysed as genes),\n",
  "         (c) blanks passed eligibility somewhere (MIN_MOL may be too high for them).")

# ---- NULL 2: position-based, for daughter_only / parent_only rows --------------------
# These have no shift; the claim is that the gene is POSITIONED inside one cell type.
# Compare |delta| against the blanks' |delta| in the same lineage.
thr_pos <- res %>% filter(is_blank, design %in% c("daughter_only", "parent_only")) %>%
  mutate(d = coalesce(delta_d9, delta_d0)) %>% filter(is.finite(d)) %>%
  group_by(lineage, design) %>%
  summarise(n_blank = n(), thr_absdelta = quantile(abs(d), .95, na.rm = TRUE), .groups = "drop")
cat("\nblank null for daughter_only / parent_only:\n"); print(as.data.frame(thr_pos), row.names = FALSE)

# ---- classify ------------------------------------------------------------------------
d2 <- res %>% filter(!is_blank, design == "D2_persist") %>%
  select(gene, lineage, d2_shift = median_shift, d2_consistency = consistency)

hits <- res %>% filter(!is_blank) %>%
  left_join(thr_d1  %>% select(-n_blank), by = "lineage") %>%
  left_join(thr_pos %>% select(-n_blank), by = c("lineage", "design")) %>%
  left_join(d2, by = c("gene", "lineage")) %>%
  mutate(
    clears_null = coalesce(case_when(
      design == "D1_lineage" ~ is.finite(consistency) & consistency >= thr_cons &
                               abs(median_shift) >= thr_shift,
      design %in% c("daughter_only", "parent_only") ~
        is.finite(thr_absdelta) & abs(coalesce(delta_d9, delta_d0)) >= thr_absdelta,
      TRUE ~ FALSE), FALSE),
    designs_agree = is.finite(d2_shift) & is.finite(median_shift) &
                    sign(d2_shift) == sign(median_shift),
    class = case_when(
      design == "parent_only"                                        ~ "lost_in_daughter",
      design == "daughter_only" &  clears_null                       ~ "daughter_positioned",
      design == "daughter_only" & !clears_null                       ~ "daughter_diffuse",
      design != "D1_lineage"                                         ~ NA_character_,
      !clears_null & abs(delta_d0) >= thr_shift &
        sign(delta_d0) == sign(delta_d9)                             ~ "preserved",
      clears_null & sign(delta_d0) == sign(delta_d9) &
        abs(delta_d9) > abs(delta_d0)                                ~ "amplified",
      clears_null & abs(delta_d0) <  thr_shift &
        abs(delta_d9) >= thr_shift                                   ~ "gained",
      clears_null & abs(delta_d0) >= thr_shift &
        abs(delta_d9) <  thr_shift                                   ~ "lost",
      clears_null                                                    ~ "shifted",
      TRUE                                                           ~ "none"),
    # Ranking: clears its null, and (for D1) both designs agree in sign.
    # `marker` is NOT a filter - it is reported as an annotation only. The clustering
    # used expression, never position, so a defining gene's POSITION inside its own type
    # is not circular. (Excluding markers would have dropped Mmp10/Mmp13/Plau - our
    # strongest candidates and the paper's own p.2017 result.)
    rank_ok = coalesce(clears_null &
                       (designs_agree | design != "D1_lineage"), FALSE)) %>%
  arrange(desc(rank_ok), desc(clears_null), desc(consistency), desc(abs(median_shift)))
write_csv(hits, file.path(OUT, "hits.csv"))

cat("\nclass x lineage (rows clearing their null):\n")
print(table(hits$class[hits$rank_ok], hits$lineage[hits$rank_ok], useNA = "ifany"))
cat(sprintf("\n%d rows | %d clear their null | %d ranked (designs agree) | %d of those flagged marker\n",
            nrow(hits), sum(hits$clears_null, na.rm = TRUE), sum(hits$rank_ok, na.rm = TRUE),
            sum(hits$rank_ok & hits$marker, na.rm = TRUE)))
cat("\ntop 20 ranked candidates:\n")
print(as.data.frame(hits %>% filter(rank_ok) %>%
  select(gene, lineage, design, class, marker, delta_d0, delta_d9, median_shift, consistency, n_pairs) %>%
  mutate(across(where(is.numeric), ~round(., 4))) %>% head(20)), row.names = FALSE)
