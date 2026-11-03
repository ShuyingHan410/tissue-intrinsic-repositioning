# Step 3 (v2): meta-analytic inference on the per-slice deltas produced by 01_run_chunk.R.
# Replaces the blank-threshold ranking in 02_merge.R. See PLAN_v2.md.
#
#   Step 2  pool slices within each arm  -> arms.csv          (REML + Knapp-Hartung)
#   Step 3  contrast the arms            -> contrasts.csv     (D1, D2, type difference + CIs)
#   Step 4  BH within lineage            -> d1_q in contrasts.csv
#   Step 5  blank z-score diagnostic     -> blank_diagnostics.csv
#
#   Rscript 03_meta.R
OUT <- Sys.getenv("OUT", ".")
suppressMessages({library(dplyr); library(readr); library(tidyr); library(purrr); library(metafor)})

RES <- file.path(OUT, "results")
slc <- list.files(RES, "^slices_.*[.]csv$", full.names = TRUE) %>%
  map_dfr(read_csv, show_col_types = FALSE, progress = FALSE)
gi  <- read_csv(file.path(OUT, "genes_index.csv"), show_col_types = FALSE)
stopifnot("se_delta" %in% names(slc))   # 01_run_chunk.R must have been re-run with v2 delta_core.R

# side -> arm label; `day` already distinguishes the two parent arms
ARM <- c(parent_D0 = "parent_D0", parent_D9 = "parent_D9", daughter_D9 = "daughter_D9")
CELL <- tribble(
  ~lineage,         ~side,          ~cell_type,
  "Fibro6->IAF2",   "parent_D0",    "Fibro 6",   "Fibro6->IAF2",   "parent_D9",   "Fibro 6",
  "Fibro6->IAF2",   "daughter_D9",  "IAF 2",     "Fibro6->IAF3",   "parent_D0",   "Fibro 6",
  "Fibro6->IAF3",   "parent_D9",    "Fibro 6",   "Fibro6->IAF3",   "daughter_D9", "IAF 3",
  "Fibro2->IAF4",   "parent_D0",    "Fibro 2",   "Fibro2->IAF4",   "parent_D9",   "Fibro 2",
  "Fibro2->IAF4",   "daughter_D9",  "IAF 4",     "Fibro7->Fibro4", "parent_D0",   "Fibro 7",
  "Fibro7->Fibro4", "parent_D9",    "Fibro 7",   "Fibro7->Fibro4", "daughter_D9", "Fibro 4")

dat <- slc %>%
  filter(eligible, is.finite(delta), is.finite(se_delta), se_delta > 0) %>%
  mutate(arm = ARM[side]) %>% filter(!is.na(arm)) %>%
  left_join(CELL, by = c("lineage", "side")) %>%
  left_join(gi %>% select(gene, is_blank), by = "gene")

# ---- STEP 2: pool slices within an arm -------------------------------------------------
# REML for tau^2 (less biased than DerSimonian-Laird), Knapp-Hartung for the interval
# (rescales by observed weighted scatter and uses t_{k-1}, so it widens most when k is small).
MIN_K <- as.integer(Sys.getenv("MIN_K", "3"))
pool <- function(y, v) {
  k <- length(y)
  if (k < MIN_K) return(NULL)
  fit <- tryCatch(rma(yi = y, vi = v, method = "REML", test = "knha"),
                  error = function(e) tryCatch(rma(yi = y, vi = v, method = "DL", test = "knha"),
                                               error = function(e2) NULL))
  if (is.null(fit)) return(NULL)
  # KH inflation factor q, reported so a reader can see whether slices scattered more or less
  # than their weights predicted
  w <- 1 / (v + fit$tau2); mu <- sum(w * y) / sum(w)
  tibble(k = k, est = as.numeric(fit$beta), se = fit$se,
         ci_lo = fit$ci.lb, ci_hi = fit$ci.ub,
         tau = sqrt(fit$tau2), I2 = fit$I2,
         q_kh = sum(w * (y - mu)^2) / (k - 1),
         p = fit$pval)
}

arms <- dat %>% group_by(gene, lineage, arm, cell_type, is_blank) %>%
  group_modify(~ { r <- pool(.x$delta, .x$se_delta^2); if (is.null(r)) tibble() else r }) %>%
  ungroup()
write_csv(arms, file.path(OUT, "arms.csv"))
cat(sprintf("arms.csv: %s rows (%d genes x lineage x arm)\n",
            format(nrow(arms), big.mark = ","), nrow(arms)))

# ---- STEP 3: contrasts, Welch-Satterthwaite df -----------------------------------------
# Arms are disjoint sets of slices, so their estimates are independent and SEs add in
# quadrature. Each arm carries its own t df (k-1), hence Welch-Satterthwaite for the pair.
contrast <- function(a, b) {
  if (nrow(a) != 1 || nrow(b) != 1) return(tibble(est = NA_real_))
  est <- a$est - b$est; se <- sqrt(a$se^2 + b$se^2)
  df  <- (a$se^2 + b$se^2)^2 / (a$se^4 / (a$k - 1) + b$se^4 / (b$k - 1))
  tc  <- qt(.975, df)
  tibble(est = est, se = se, ci_lo = est - tc * se, ci_hi = est + tc * se,
         df = df, p = 2 * pt(-abs(est / se), df))
}
grab <- function(d, a) d %>% filter(arm == a)

contrasts <- arms %>% group_by(gene, lineage, is_blank) %>% group_modify(function(d, key) {
  p0 <- grab(d, "parent_D0"); p9 <- grab(d, "parent_D9"); d9 <- grab(d, "daughter_D9")
  c1 <- contrast(d9, p0); c2 <- contrast(p9, p0); ct <- contrast(d9, p9)
  tibble(
    k_d0 = if (nrow(p0)) p0$k else NA_integer_,
    k_p9 = if (nrow(p9)) p9$k else NA_integer_,
    k_d9 = if (nrow(d9)) d9$k else NA_integer_,
    # arm-level estimates WITH intervals: these are what daughter_only / parent_only report,
    # since for those routes no cross-day contrast is meaningful
    delta_d0 = if (nrow(p0)) p0$est else NA_real_,
    d0_ci_lo = if (nrow(p0)) p0$ci_lo else NA_real_,
    d0_ci_hi = if (nrow(p0)) p0$ci_hi else NA_real_,
    d0_p     = if (nrow(p0)) p0$p else NA_real_,
    delta_d9 = if (nrow(d9)) d9$est else NA_real_,
    d9_ci_lo = if (nrow(d9)) d9$ci_lo else NA_real_,
    d9_ci_hi = if (nrow(d9)) d9$ci_hi else NA_real_,
    d9_p     = if (nrow(d9)) d9$p   else NA_real_,
    I2_daughter = if (nrow(d9)) d9$I2 else NA_real_,
    tau_daughter = if (nrow(d9)) d9$tau else NA_real_,
    d1 = c1$est, d1_se = c1$se, d1_ci_lo = c1$ci_lo, d1_ci_hi = c1$ci_hi, d1_df = c1$df, d1_p = c1$p,
    d2 = c2$est, d2_se = c2$se, d2_ci_lo = c2$ci_lo, d2_ci_hi = c2$ci_hi, d2_df = c2$df, d2_p = c2$p,
    typediff = ct$est, typediff_se = ct$se,
    typediff_ci_lo = ct$ci_lo, typediff_ci_hi = ct$ci_hi, typediff_p = ct$p)
}) %>% ungroup()

PARENT <- c("Fibro6->IAF2"="Fibro 6","Fibro6->IAF3"="Fibro 6",
            "Fibro2->IAF4"="Fibro 2","Fibro7->Fibro4"="Fibro 7")
DAUGHT <- c("Fibro6->IAF2"="IAF 2","Fibro6->IAF3"="IAF 3",
            "Fibro2->IAF4"="IAF 4","Fibro7->Fibro4"="Fibro 4")

# ROUTE comes from v1's inclusion_criteria.R (the paper's 2.813-sigma expression test), NOT
# from whether an arm happened to pool. Re-deriving it from data availability silently
# reclassified every daughter_only/parent_only gene as D1_full: a gene can be poolable in both
# arms while still being *unexpressed* at Day 0, and a D1 contrast is then meaningless.
v1route <- read_csv(file.path(OUT, "gene_spatial_lineage_results.csv"), show_col_types = FALSE) %>%
  filter(!is.na(route)) %>% distinct(gene, lineage, route) %>%
  group_by(gene, lineage) %>% slice(1) %>% ungroup()

contrasts <- contrasts %>%
  left_join(v1route, by = c("gene", "lineage")) %>%
  mutate(
    parent = PARENT[lineage], daughter = DAUGHT[lineage],
    route = coalesce(route, "none"),
    # a D1 contrast is only meaningful on the D1_full route
    across(c(d1, d1_se, d1_ci_lo, d1_ci_hi, d1_df, d1_p,
             d2, d2_se, d2_ci_lo, d2_ci_hi, d2_df, d2_p),
           ~ ifelse(route == "D1_full", .x, NA_real_)),
    direction = case_when(is.finite(d1) & d1 < 0        ~ "luminal",
                          is.finite(d1) & d1 > 0        ~ "crypt-base",
                          is.finite(delta_d9) & delta_d9 < 0 ~ "luminal",
                          is.finite(delta_d9)           ~ "crypt-base", TRUE ~ NA_character_),
    # PATTERN: read BOTH endpoints, not just their difference. The same D1 arises from very
    # different biology - a gradient appearing, intensifying, reversing or dissolving - and
    # each endpoint's own CI says whether the gene was positioned at that timepoint.
    d0_positioned = is.finite(d0_ci_lo) & !(d0_ci_lo < 0 & d0_ci_hi > 0),
    d9_positioned = is.finite(d9_ci_lo) & !(d9_ci_lo < 0 & d9_ci_hi > 0),
    pattern = case_when(
      route != "D1_full"                                            ~ NA_character_,
      !d0_positioned &  d9_positioned                               ~ "gained",
       d0_positioned & !d9_positioned                               ~ "lost",
       d0_positioned &  d9_positioned & sign(delta_d0) == sign(delta_d9) &
         abs(delta_d9) > abs(delta_d0)                              ~ "amplified",
       d0_positioned &  d9_positioned & sign(delta_d0) == sign(delta_d9) ~ "dampened",
       d0_positioned &  d9_positioned                               ~ "reversed",
      TRUE                                                          ~ "weak-both-ends"),
    # POWER GUARD for `gained`: a Day-0 CI that spans zero may mean "no gradient" OR simply
    # "too few molecules to detect one". Only the former is a real gain. Flag the ambiguous
    # case: is the D0 interval even narrow enough to have excluded a gradient of size |delta_d9|?
    d0_conclusive = ifelse(pattern == "gained",
                           pmax(abs(d0_ci_lo), abs(d0_ci_hi)) < abs(delta_d9), NA)) %>%
  ungroup()

# ---- STEP 4a: recalibrate against the blank null, THEN BH ------------------------------
# The exact SEs fix expression-dependence, but slices within a mouse remain correlated and
# tau^2 only absorbs part of that. The blanks measure the residual inflation directly: with a
# correctly-calibrated test their p-values would be uniform. They are not, so rescale the test
# statistics by lambda (genomic-control style) before any FDR.
#   lambda = median(t^2 over blanks) / median(chi^2_1) ,  median(chi^2_1) = 0.4549364
gc_lambda <- function(p, df) {
  p <- p[is.finite(p) & p > 0 & p < 1]; if (length(p) < 30) return(1)
  max(1, median(qchisq(p, df = 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(.5, 1))
}
lam_d1 <- gc_lambda(contrasts$d1_p[contrasts$is_blank])
lam_d9 <- gc_lambda(contrasts$d9_p[contrasts$is_blank])
lam_d0 <- gc_lambda(contrasts$d0_p[contrasts$is_blank])
cat(sprintf("\nblank recalibration (Step 4a): lambda_D1 = %.3f  lambda_D9arm = %.3f  lambda_D0arm = %.3f\n",
            lam_d1, lam_d9, lam_d0))
recal <- function(p, lam) {
  ifelse(is.finite(p), pchisq(qchisq(p, 1, lower.tail = FALSE) / lam, 1, lower.tail = FALSE), NA_real_)
}
contrasts <- contrasts %>%
  mutate(d1_p_raw = d1_p, d9_p_raw = d9_p, d0_p_raw = d0_p,
         d1_p = recal(d1_p, lam_d1), d9_p = recal(d9_p, lam_d9), d0_p = recal(d0_p, lam_d0),
         # ONE headline test per gene, chosen by route: D1_full asks "did it move?",
         # daughter_only / parent_only ask "is it positioned?" in the arm where it exists.
         effect    = case_when(route == "D1_full" ~ d1,
                               route == "daughter_only" ~ delta_d9,
                               route == "parent_only"   ~ delta_d0, TRUE ~ NA_real_),
         effect_lo = case_when(route == "D1_full" ~ d1_ci_lo,
                               route == "daughter_only" ~ d9_ci_lo,
                               route == "parent_only"   ~ d0_ci_lo, TRUE ~ NA_real_),
         effect_hi = case_when(route == "D1_full" ~ d1_ci_hi,
                               route == "daughter_only" ~ d9_ci_hi,
                               route == "parent_only"   ~ d0_ci_hi, TRUE ~ NA_real_),
         p         = case_when(route == "D1_full" ~ d1_p,
                               route == "daughter_only" ~ d9_p,
                               route == "parent_only"   ~ d0_p, TRUE ~ NA_real_)) %>%
  # STEP 4b: BH within lineage x route. Routes ask different questions and have very different
  # pool sizes (328 D1_full vs 1350 parent_only), so they are not one exchangeable family.
  group_by(lineage, route) %>%
  mutate(q = ifelse(is_blank, NA_real_, p.adjust(ifelse(is_blank, NA_real_, p), "BH"))) %>%
  ungroup() %>%
  # expression context: without it a reader cannot tell whether a shift accompanied a large
  # abundance change. delta itself is scale-invariant (Cov(n,r)/mean(n)), and empirically
  # |D1| is uncorrelated with log2FC, but the numbers belong in the table.
  #
  # NORMALISATION MATTERS. Molecules-per-1000-CELLS is NOT comparable across these arms: the
  # daughter populations capture ~40% fewer transcripts per cell than their parents
  # (IAF3/Fibro6 = 0.60x, IAF4/Fibro2 = 0.60x, Fibro4/Fibro7 = 0.69x), so a per-cell rate
  # depresses every daughter arm by ~0.7 in log2 and 90 of 328 genes flip sign. The primary
  # column is therefore per-10k-TRANSCRIPTS (depth-normalised); the per-cell version is kept
  # alongside it, since which one is right depends on the question (per-cell = absolute output,
  # per-transcript = share of the transcriptome).
  left_join(
    slc %>% filter(eligible, side %in% c("parent_D0", "daughter_D9")) %>%
      group_by(gene, lineage, side) %>%
      summarise(rate_cell = median(rate_per_1k_cells),
                rate_tx   = median(rate_per_10k_tx), .groups = "drop") %>%
      pivot_wider(names_from = side, values_from = c(rate_cell, rate_tx)) %>%
      rename(rate_cell_d0 = rate_cell_parent_D0, rate_cell_d9 = rate_cell_daughter_D9,
             rate_tx_d0   = rate_tx_parent_D0,   rate_tx_d9   = rate_tx_daughter_D9),
    by = c("gene", "lineage")) %>%
  mutate(log2fc_expr      = log2((rate_tx_d9   + 0.1) / (rate_tx_d0   + 0.1)),
         log2fc_percell   = log2((rate_cell_d9 + 0.1) / (rate_cell_d0 + 0.1))) %>%
  select(gene, lineage, parent, daughter, route, pattern, direction,
         effect, effect_lo, effect_hi, p, q,                 # headline, route-appropriate
         delta_d0, d0_ci_lo, d0_ci_hi, d0_p, d0_positioned,
         delta_d9, d9_ci_lo, d9_ci_hi, d9_p, d9_positioned, d0_conclusive,
         rate_tx_d0, rate_tx_d9, log2fc_expr, rate_cell_d0, rate_cell_d9, log2fc_percell,
         k_d0, k_p9, k_d9, I2_daughter, tau_daughter,
         d1, d1_se, d1_ci_lo, d1_ci_hi, d1_df, d1_p, d1_p_raw,
         # D2 and the type difference are RETAINED but demoted: the surviving parent
         # population at Day 9 sits in systematically less-ulcerated tissue (median 4.9% vs
         # 0.8% epithelium remaining, p = 0.016; absent from 12 of 20 daughter slices), so it
         # is NOT a valid control for the daughter arm. Descriptive only - do not grade D1 on it.
         d2, d2_se, d2_ci_lo, d2_ci_hi, d2_df, d2_p,
         starts_with("typediff"), is_blank) %>%
  arrange(q, desc(abs(effect)))

write_csv(contrasts %>% filter(!is_blank) %>% select(-is_blank),
          file.path(OUT, "contrasts.csv"))
cat(sprintf("contrasts.csv: %s gene x lineage rows\n",
            format(sum(!contrasts$is_blank), big.mark = ",")))

# ---- STEP 5: blank diagnostics ---------------------------------------------------------
# The SEs above are only as good as the Step-1 formula. Blanks carry no transcript, so their
# z MUST be ~N(0,1). If not, everything above is mis-scaled. THIS GATES THE RESULT.
bd <- slc %>% left_join(gi %>% select(gene, is_blank), by = "gene") %>%
  filter(is_blank, eligible, is.finite(z_delta))
write_csv(bd %>% select(gene, lineage, side, Slice_ID, cells, mol, delta, se_delta, z_delta),
          file.path(OUT, "blank_diagnostics.csv"))
lam <- stats::sd(bd$z_delta, na.rm = TRUE)
cat(sprintf("\nblank_diagnostics.csv: %s z-scores from %d blank barcodes\n",
            format(nrow(bd), big.mark = ","), n_distinct(bd$gene)))
cat(sprintf("  mean(z) = %+.3f   (target 0)\n", mean(bd$z_delta, na.rm = TRUE)))
cat(sprintf("  sd(z)   = %.3f    (target 1)\n", lam))
cat(sprintf("  %%|z|>1.96 = %.1f%%  (target 5.0%%)\n", 100 * mean(abs(bd$z_delta) > 1.96, na.rm = TRUE)))
if (lam > 1.15 || lam < 0.85)
  cat(sprintf("  *** SE MIS-CALIBRATED: divide all z by lambda = %.3f before interpreting ***\n", lam))
if (nrow(bd) < 100)
  cat("  *** too few blank z-scores to judge calibration ***\n")

# blank q-values should be uniform; a spike near 0 means the model is mis-specified
bq <- contrasts %>% filter(is_blank, is.finite(d1_p))
if (nrow(bq)) cat(sprintf("  blank D1 p-values below 0.05: %.1f%% (target ~5%%, n=%d)\n",
                          100 * mean(bq$d1_p < 0.05), nrow(bq)))

# ---- summary ---------------------------------------------------------------------------
cc <- contrasts %>% filter(!is_blank)
cat("\nroute x pattern:\n"); print(table(cc$route, cc$pattern, useNA = "ifany"))
cat("\nsignificant (q < 0.05) per route:\n")
print(cc %>% filter(q < 0.05) %>% count(route, lineage) %>% as.data.frame(), row.names = FALSE)
cat("\ntop 20 by q:\n")
print(as.data.frame(cc %>% filter(q < 0.05) %>%
  transmute(gene, lineage, route, pattern,
            effect = round(effect, 4),
            ci = sprintf("[%+.4f, %+.4f]", effect_lo, effect_hi),
            q = signif(q, 2), I2 = round(I2_daughter),
            log2FC = round(log2fc_expr, 2)) %>%
  head(20)), row.names = FALSE)
ng <- sum(cc$route == "D1_full" & cc$q < 0.05 & cc$pattern == "gained", na.rm = TRUE)
nc <- sum(cc$route == "D1_full" & cc$q < 0.05 & cc$pattern == "gained" & !cc$d0_conclusive, na.rm = TRUE)
cat(sprintf("\n'gained' calls where the Day-0 CI is too wide to rule out a gradient: %d of %d\n", nc, ng))
