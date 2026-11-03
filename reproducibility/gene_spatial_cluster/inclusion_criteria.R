# =====================================================================================
# INCLUSION / EXCLUSION CRITERIA for the gene-spatial-lineage sweep.
# Sourced by delta_core.R. Set these BEFORE looking at results - they are fixed by the
# blanks and the coverage counts, never tuned on the hit list.
#
# Every gene x lineage gets a ROUTE, never a silent drop:
#   D1_full        expressed in parent@D0 AND daughter@D9 -> parent->daughter shift (D1+D2)
#   daughter_only  absent/sparse at D0, expressed at D9   -> within-daughter delta vs blanks
#   parent_only    expressed at D0, lost at D9            -> within-parent delta + level change
#   excluded       hard-excluded, or not expressed either side
# =====================================================================================

# ---- 1. HARD EXCLUSIONS -------------------------------------------------------------
# The authors removed these as false-positive contaminated (STAR Methods): "while they
# were included in clustering, the expression of some of these genes may be contaminated
# by false positives, and we removed them from subsequent analysis."  All six are still
# present in adata.h5ad, so we must drop them ourselves.
PAPER_FP_REMOVED <- c("Irf8", "Foxl1", "Foxf1", "Cyba", "Foxp3", "Il17a")

# Blank barcodes are never results - they are the null. Detected by prefix.
#' Blank status comes from genes_index.csv, NOT from the barcode's name.
#' The h5ad carries an earlier panel annotation in which two control barcodes are named
#' `Dclk1` and `Gdf1`; the published Table S2 lists those same positions as Blank-48 and
#' Blank-49. A `^Blank` name match therefore found only 48 of the 50 controls, shrinking the
#' null and analysing two control probes as if they were genes. genes_index.csv is now
#' corrected to the paper's 940 genes + 50 blanks, and is the single source of truth.
BLANK_SET <- local({
  f <- file.path(Sys.getenv("OUT", "."), "genes_index.csv")
  gi <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  s <- gi$gene[tolower(as.character(gi$is_blank)) %in% c("true", "1")]
  if (length(s) != 50) warning(sprintf(
    "expected 50 blank barcodes, genes_index.csv has %d - the blank null will be mis-sized", length(s)))
  s
})
is_blank_gene <- function(g) g %in% BLANK_SET

# ---- 2. MARKER ANNOTATION (reporting only - NOT a filter) ----------------------------
# Genes the paper names as defining each lineage's types (Fig 7A, p.2013 marker lists).
# Recorded as a column so a reader can see it, but it does NOT gate testing or ranking:
# the clustering used EXPRESSION and never position, so a defining gene's POSITION inside
# its own type is not circular. Every gene passes through the same criteria.
LINEAGE_SIGNATURE <- list(
  "Fibro6->IAF2"   = c("Igfbp5","Grem1","Col18a1","Sdc1","Col1a2","C3","Vcam1",
                       "Wnt2b","Col15a1","Cd34"),
  "Fibro6->IAF3"   = c("Igfbp5","Grem1","Col18a1","Mmp2","Sdc1","Col1a2","C3","Vcam1",
                       "Wnt2b","Col15a1","Cd34"),
  "Fibro2->IAF4"   = c("Il11","Il1b","Il1rl1","Mmp3","Mmp10","Mmp13","Plau",
                       "Adamdec1","Wnt5a","Bmp2","Bmp5","Bmp7","Pdgfra"),
  "Fibro7->Fibro4" = c("Il1b","Il11","Csf3","Cxcl5","Mmp3","Pi16","Ackr4","Col14a1",
                       "Grem1","Cd34","C3","Cd248")
)
is_defining_marker <- function(gene, lineage) {
  sig <- LINEAGE_SIGNATURE[[lineage]]
  !is.null(sig) && gene %in% sig
}

# ---- 3. EXPRESSION TEST (the paper's own criterion) ----------------------------------
# Cadinu et al. classify a gene as 'expressed' in a cluster using the 48/50 blank
# barcodes within that cluster and the one-sided Vysochanskij-Petunin inequality:
#   "a gene must have an expression value that is 2.813s greater than m in order to
#    meet threshold" (5% probability cutoff).
# We apply the identical rule per CELL TYPE, on per-slice totals.
VP_K         <- 2.813
MIN_COVERAGE <- as.numeric(Sys.getenv("MIN_COVERAGE", "0.5"))  # >=50% of available slices
MIN_SLICES   <- as.integer(Sys.getenv("MIN_SLICES",   "3"))    # >=3 eligible SLICES per side
MIN_RATE_PER_1K <- as.numeric(Sys.getenv("MIN_RATE_PER_1K", "8"))  # retained for reporting only (v16 gate: SE0)
if (!exists("SE_MAX")) SE_MAX <- as.numeric(Sys.getenv("SE_MAX", "0.05"))  # per-slice precision gate, frozen 2026-09-13

#' dat  : joined table (Slice_ID, Tier2, day, mouse, r, <gene cols>, blank_counts)
#' Returns TRUE if the gene's per-slice total in this type x day exceeds
#'   mean(blank barcodes) + 2.813 * sd(blank barcodes)
#' in at least MIN_COVERAGE of slices.
#'
#' FIXED 2026-08-17. The previous version computed
#'     thr <- mean(per$b_mean) + VP_K * sd(per$b_mean)
#' where b_mean was one value PER SLICE (that slice's blank total / 48). That takes the
#' SD *across slices*, which measures how much slices differ in size and depth - not the
#' spread across the 48 blank barcodes that the Vysochanskij-Petunin rule calls for.
#' Slices vary ~100-fold, so the across-slice SD swamped the mean (IAF 4 @ D9: mean 28.8,
#' sd 32.2 -> thr 119.4) and genuinely expressed genes were routed to `excluded`:
#' Il1b (13.7x above a single blank barcode, 18/18 slices), Cxcl5, and Il11 in the SM
#' lineage. Per-barcode thresholds are now precomputed by 00b_blank_thresholds.py, which
#' has access to all 48 blank columns (a chunk job sees only its own 30 genes).
BLANK_THR <- local({
  f <- file.path(Sys.getenv("OUT", "."), "blank_thresholds.csv")
  if (!file.exists(f)) stop(
    "blank_thresholds.csv not found - run `python3 00b_blank_thresholds.py` first.\n",
    "  It supplies the per-(slice x cell type x day) 2.813-sigma detection threshold\n",
    "  computed across the 48 blank barcodes, which the chunk jobs cannot derive alone.")
  suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
})

#' EXPR_SCOPE = "cluster" (default) reproduces the paper exactly: one bar per cluster,
#'   from each blank barcode's average count per cell pooled over every cell of that type,
#'   and the gene's own average count per cell compared against it. No slice-consistency
#'   requirement - the paper imposes none.
#' EXPR_SCOPE = "slice" is the stricter local variant: each slice judged against its own
#'   bar, gene must clear it in >=MIN_COVERAGE of slices. Guards against one slice carrying
#'   a gene, but is not what the paper describes and runs 3-37% harsher.
#' EXPR_SCOPE = "tableS2" (DEFAULT) uses the paper's PUBLISHED per-cluster gene lists
#'   (Table S2, mmc2.xlsx -> tableS2_expressed.csv). No threshold is recomputed: a gene counts
#'   as expressed in a cell type iff the authors list it. This is exactly reproducible and
#'   sidesteps the raw-vs-normalised discrepancy - their expression values are scanpy
#'   normalize_total + log10, ours are raw counts per cell, which is why our own
#'   reimplementation passed 18-49 extra genes per cluster.
#'   CAVEAT: Table S2 pools all timepoints per cluster, so it cannot say whether a gene is
#'   present at Day 0 vs Day 9. For the DAUGHTER types that is harmless - IAF 2/3/4 are
#'   essentially Day-9-only populations (IAF 4 has 10 cells at Day 0). For the PARENT types,
#'   which exist at both days, the list is dominated by the more abundant Day-0 cells. We
#'   therefore AND it with the per-day test so the route still distinguishes the timepoints.
EXPR_SCOPE <- Sys.getenv("EXPR_SCOPE", "tableS2")

TABLE_S2 <- local({
  f <- file.path(Sys.getenv("OUT", "."), "tableS2_expressed.csv")
  if (EXPR_SCOPE == "tableS2" && !file.exists(f)) stop(
    "tableS2_expressed.csv not found - export it from the paper's mmc2.xlsx first.\n",
    "  It is the authors' published list of genes passing their 2.813-sigma bar per cluster.")
  if (file.exists(f)) suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  else NULL
})

expressed_in <- function(dat, type, gene, dy, ...) {
  z <- dat %>% filter(Tier2 == type, day == dy)
  if (!nrow(z)) return(NA)

  if (EXPR_SCOPE == "tableS2") {
    if (!type %in% TABLE_S2$Tier2) return(NA)          # cluster absent from Table S2
    in_s2 <- any(TABLE_S2$Tier2 == type & TABLE_S2$gene == gene)
    if (!in_s2) return(FALSE)
    # AND the per-day cluster test, so Day 0 vs Day 9 is still distinguished
    b <- BLANK_THR %>% filter(scope == "cluster", Tier2 == type, day == dy)
    if (!nrow(b) || b$cells[1] < MIN_CELLS) return(NA)
    return(mean(z[[gene]]) >= b$thr_per_cell[1])
  }

  if (EXPR_SCOPE == "cluster") {
    b <- BLANK_THR %>% filter(scope == "cluster", Tier2 == type, day == dy)
    if (!nrow(b) || b$cells[1] < MIN_CELLS) return(NA)
    # gene's average count per cell across all cells of this type x day
    return(mean(z[[gene]]) >= b$thr_per_cell[1])
  }

  per <- z %>% group_by(Slice_ID) %>%
    summarise(cells = n(), g = sum(.data[[gene]]), .groups = "drop") %>%
    filter(cells >= MIN_CELLS) %>%
    inner_join(BLANK_THR %>% filter(scope == "slice", Tier2 == type, day == dy) %>%
                 select(Slice_ID, thr), by = "Slice_ID")
  if (nrow(per) < 2) return(NA)
  mean(per$g >= per$thr, na.rm = TRUE) >= MIN_COVERAGE
}

# ---- 4. COVERAGE + REPLICATION (constants declared above, used here) -----------------
# NOTE: SLICE-PAIR design - every eligible Day0 slice x every eligible Day9 slice
# (up to 10 x 21 = 210 in MU). Mouse IDs are carried in the pairs table so the
# dependence can be modelled or ignored downstream, as you choose.

#' Slice-level coverage for one gene in one type x day: n_eligible / n_available
coverage_of <- function(dat, type, gene, dy) {
  z <- dat %>% filter(Tier2 == type, day == dy) %>%
    group_by(Slice_ID) %>%
    summarise(cells = n(), mol = sum(.data[[gene]]), sd_r = stats::sd(r),
              sum_w2 = ifelse(sum(.data[[gene]]) > 0,
                              sum(.data[[gene]]^2) / sum(.data[[gene]])^2 - 1 / n(), NA_real_),
              .groups = "drop") %>%
    mutate(se0 = ifelse(is.finite(sum_w2) & sum_w2 > 0 & is.finite(sd_r), sd_r * sqrt(sum_w2), NA_real_))
  avail <- sum(z$cells >= MIN_CELLS)
  # same per-slice gate as delta_core.R (frozen 2026-09-13, plan v16 §0.6): cells AND SE0 < SE_MAX
  elig  <- sum(z$cells >= MIN_CELLS & is.finite(z$se0) & z$se0 > 0 & z$se0 < SE_MAX)
  list(n_avail = avail, n_elig = elig, frac = if (avail) elig / avail else 0,
       med_mol = if (avail) stats::median(z$mol[z$cells >= MIN_CELLS]) else NA_real_)
}

# ---- 5. ROUTING ----------------------------------------------------------------------
#' Decide, for one gene x lineage, which analysis it can support.
#' Returns list(route, exclude_reason, marker, plus coverage on both sides).
route_gene <- function(dat, gene, lg) {
  if (gene %in% PAPER_FP_REMOVED)
    return(list(route = "excluded", exclude_reason = "paper_FP_removed", marker = NA,
                p_elig = NA, p_avail = NA, p_med_mol = NA,
                d_elig = NA, d_avail = NA, d_med_mol = NA))
  # BLANKS ARE THE NULL, NOT CANDIDATES. They must go through the SAME delta machinery as
  # real genes, but they cannot be asked to pass the expression test they themselves
  # define (a blank cannot exceed the blank background), nor a molecule floor written for
  # real genes (a blank carries ~0.034 counts/cell). Route them on coverage alone.
  if (is_blank_gene(gene)) {
    cp <- coverage_of(dat, lg$parent,   gene, "Day 0")
    cd <- coverage_of(dat, lg$daughter, gene, "Day 9")
    ok_p <- cp$n_elig >= MIN_SLICES; ok_d <- cd$n_elig >= MIN_SLICES
    route <- if (ok_p && ok_d) "D1_full" else if (!ok_p && ok_d) "daughter_only" else
             if (ok_p && !ok_d) "parent_only" else "excluded"
    return(list(route = route, exclude_reason = if (route=="excluded") "blank_too_sparse" else NA_character_,
                marker = FALSE,
                p_elig = cp$n_elig, p_avail = cp$n_avail, p_med_mol = cp$med_mol,
                d_elig = cd$n_elig, d_avail = cd$n_avail, d_med_mol = cd$med_mol))
  }
  cp <- coverage_of(dat, lg$parent,   gene, "Day 0")
  cd <- coverage_of(dat, lg$daughter, gene, "Day 9")
  ok_p <- cp$frac >= MIN_COVERAGE && isTRUE(expressed_in(dat, lg$parent,   gene, "Day 0"))
  ok_d <- cd$frac >= MIN_COVERAGE && isTRUE(expressed_in(dat, lg$daughter, gene, "Day 9"))
  route <- if (ok_p && ok_d) "D1_full" else if (!ok_p && ok_d) "daughter_only" else
           if (ok_p && !ok_d) "parent_only" else "excluded"
  list(route = route,
       exclude_reason = if (route == "excluded") "not_expressed_either_side" else NA_character_,
       marker = is_defining_marker(gene, lg$name),      # ranking flag, NOT an exclusion
       p_elig = cp$n_elig, p_avail = cp$n_avail, p_med_mol = cp$med_mol,
       d_elig = cd$n_elig, d_avail = cd$n_avail, d_med_mol = cd$med_mol)
}
