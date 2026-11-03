# Shared core for the gene-spatial-lineage sweep. Sourced by 01_run_chunk.R (cluster)
# and test_local.R (laptop). No cluster-specific code in here.
#
# The statistic, per gene x cell type x slice   (validated: Wnt2b/Bmp2 opposite signs;
# within-IAF4 Mmp13/Mmp10/Plau luminal, p = 6e-7; blanks ~56% consistency):
#     delta = sum(n * r)/sum(n) - mean(r over ALL cells of the type in that slice)
# Orientation of r checked per slice via TA vs Colonocytes (low r = crypt base).
suppressMessages({library(dplyr); library(readr); library(tidyr)})

MIN_CELLS <- as.integer(Sys.getenv("MIN_CELLS", "30"))
MIN_MOL   <- as.integer(Sys.getenv("MIN_MOL",   "20"))   # retained for reporting only; no longer gates
SE_MAX    <- as.numeric(Sys.getenv("SE_MAX",   "0.05"))  # per-slice precision gate (v16 §0.6, frozen 2026-09-13)

# inclusion/exclusion criteria: hard exclusions, circularity screen, the paper's
# 2.813-sigma expression test, coverage, and per-gene ROUTING (D1_full /
# daughter_only / parent_only / excluded).  Fixed before results are inspected.
source(file.path(Sys.getenv("OUT", "."), "inclusion_criteria.R"))

# Design 1 = parent(D0) vs daughter(D9): across the lineage (composition-proof).
# Design 2 = parent(D0) vs parent(D9): same label persisting (selection-proof).
LINEAGES <- list(
  list(name = "Fibro6->IAF2",   parent = "Fibro 6", daughter = "IAF 2",   layer = "MU"),
  list(name = "Fibro6->IAF3",   parent = "Fibro 6", daughter = "IAF 3",   layer = "MU"),
  list(name = "Fibro2->IAF4",   parent = "Fibro 2", daughter = "IAF 4",   layer = "MU"),
  list(name = "Fibro7->Fibro4", parent = "Fibro 7", daughter = "Fibro 4", layer = "SM")
)

# ---- coordinates: one layer's *_xyt_mapped.csv files -> Slice_ID, Cell_ID, r ----------
load_coords <- function(coord_dir, layer = "MU") {
  fs <- list.files(coord_dir, sprintf("_%s_xyt_mapped[.]csv$", layer), full.names = TRUE)
  if (!length(fs)) return(NULL)
  bind_rows(lapply(fs, function(f) {
    d <- suppressMessages(read_csv(f, show_col_types = FALSE, progress = FALSE))
    d %>% transmute(Slice_ID = as.character(Slice_ID), Cell_ID = as.character(Cell_ID),
                    r = as.numeric(r_for_jsd))
  })) %>% filter(is.finite(r), r >= 0, r <= 1)
}

# ---- per-slice orientation check, HARMONISED convention (2026-08-10): inner(lumen)=0 ->
#      outer=1, so Colonocytes (luminal surface) must sit BELOW TA (crypt base = deep). ----
# Emits ONE ROW PER SLICE IN `coords`, never fewer. Building the table from the epithelial
# cells instead made slices with no surviving epithelium vanish silently, so the row count
# read as a slice count and under-reported coverage (27 rows for 31 analysed slices).
# Ulcerated Day-9 sections legitimately have no TA/Colonocytes left inside the MU mask; that
# is "not evaluable", which is a different statement from "inverted" and must not look like one.
MIN_ORIENT_CELLS <- as.integer(Sys.getenv("MIN_ORIENT_CELLS", "30"))
orientation_table <- function(meta, coords) {
  ep <- meta %>% filter(Tier2 %in% c("TA", "Colonocytes")) %>%
    inner_join(coords, by = c("Slice_ID", "Cell_ID")) %>%
    group_by(Slice_ID, Tier2) %>% summarise(n = n(), m = mean(r), .groups = "drop") %>%
    # a mean over a handful of cells is noise, not a verdict
    mutate(m = ifelse(n >= MIN_ORIENT_CELLS, m, NA_real_)) %>%
    pivot_wider(names_from = Tier2, values_from = c(n, m), values_fill = list(n = 0L))
  for (nm in c("n_TA", "n_Colonocytes", "m_TA", "m_Colonocytes"))
    if (!nm %in% names(ep)) ep[[nm]] <- if (startsWith(nm, "n_")) 0L else NA_real_
  tibble(Slice_ID = sort(unique(coords$Slice_ID))) %>%
    left_join(ep, by = "Slice_ID") %>%
    mutate(n_TA = coalesce(n_TA, 0L), n_Colonocytes = coalesce(n_Colonocytes, 0L),
           TA = m_TA, Colonocytes = m_Colonocytes,
           evaluable = is.finite(TA) & is.finite(Colonocytes),
           # TRUE / FALSE only when the comparison could actually be made; NA otherwise, so
           # "could not check" never masquerades as "check failed".
           orientation_ok = ifelse(evaluable, Colonocytes < TA, NA),
           status = case_when(!evaluable    ~ "not_evaluable",
                              orientation_ok ~ "ok",
                              TRUE           ~ "INVERTED")) %>%
    select(Slice_ID, n_TA, n_Colonocytes, TA, Colonocytes, evaluable, orientation_ok, status)
}

# ---- delta per slice for one gene within one type -------------------------------------
#' Per-slice delta PLUS every quantity it is built from, so the statistic can be
#' recomputed or replaced downstream without rerunning the pipeline:
#'   rbar_gene = sum(n*r)/sum(n)      rbar_cells = mean(r)      delta = rbar_gene - rbar_cells
#' also: n_pos (cells with >0 molecules), sd_r, and the molecule-weighted sd.
#'
#' v2 (Step 1): every delta now carries an EXACT standard error. Writing w_i = n_i/sum(n),
#'   delta = sum_i (w_i - 1/N) r_i
#' is a fixed linear combination of the positions, so under the null "expression is unrelated
#' to position" (i.e. permuting r among the N cells) its variance is exact, no simulation:
#'   Var0(delta) = Var(r) * N/(N-1) * sum_i (w_i - 1/N)^2
#'               = var(r) * ( sum(n^2)/sum(n)^2 - 1/N )        [var() = sample variance]
#' the second form is what is computed here (algebraically identical, verified numerically).
#' This is SELF-CALIBRATING: a sparse gene automatically gets a wider SE, which is why the
#' blank barcodes are no longer needed to set a per-lineage threshold.
delta_by_slice <- function(dat, type, gene) {
  dat %>% filter(Tier2 == type) %>%
    group_by(Slice_ID, day, mouse) %>%
    summarise(cells = n(), n_pos = sum(.data[[gene]] > 0), mol = sum(.data[[gene]]),
              sum_nr = sum(.data[[gene]] * r), sum_r = sum(r),
              rbar_gene  = ifelse(sum(.data[[gene]]) > 0, sum(.data[[gene]] * r) / sum(.data[[gene]]), NA_real_),
              rbar_cells = mean(r),
              sd_r = stats::sd(r),
              sd_r_gene = ifelse(sum(.data[[gene]]) > 1,
                sqrt(sum(.data[[gene]] * (r - sum(.data[[gene]] * r)/sum(.data[[gene]]))^2) /
                     (sum(.data[[gene]]) - 1)), NA_real_),
              # sum_i (w_i - 1/N)^2, kept so the SE can be recomputed or reweighted downstream
              sum_w2 = ifelse(sum(.data[[gene]]) > 0,
                              sum(.data[[gene]]^2) / sum(.data[[gene]])^2 - 1 / n(), NA_real_),
              tot_counts = sum(total_counts),
              rate_per_1k_cells = 1000 * sum(.data[[gene]]) / n(),
              rate_per_10k_tx   = 1e4  * sum(.data[[gene]]) / sum(total_counts),
              .groups = "drop") %>%
    # sum_w2 can be <=0 only when a single cell carries every molecule; guard so sqrt() is safe
    mutate(se_delta = ifelse(is.finite(sum_w2) & sum_w2 > 0 & is.finite(sd_r),
                             sd_r * sqrt(sum_w2), NA_real_),
           # ELIGIBILITY (frozen 2026-09-13, simulation plan v16 §0.6, option (a)):
           #   >= MIN_CELLS cells AND 0 < SE0 < SE_MAX.
           # SE0 measures the positional information the estimator actually has, so a sparse
           # or concentrated slice fails automatically. Replaces the molecule/rate floor
           # (MIN_MOL | MIN_RATE_PER_1K), which admitted 1-2-molecule slices on Day 9 and
           # was never the rule the simulation benchmark validated. Same rule as
           # 01c_stratified_chunk.R (SE_MAX) and simulation/pilot/engine.R.
           eligible = cells >= MIN_CELLS & is.finite(se_delta) & se_delta > 0 & se_delta < SE_MAX,
           delta    = ifelse(eligible, rbar_gene - rbar_cells, NA_real_),
           z_delta  = ifelse(is.finite(se_delta) & se_delta > 0, delta / se_delta, NA_real_))
}
#' the eligible subset, used for pairing
eligible_slices <- function(x) x %>% filter(eligible, is.finite(delta))

# ---- mouse pairs -> median shift + direction consistency ------------------------------
#' SLICE-PAIR design: every eligible Day0 slice is paired with every eligible Day9 slice.
#' No mouse-level collapsing - mouse IDs are carried through as columns so any grouping
#' can still be applied downstream. Returns list(summary, pairs).
pair_stats <- function(p_d0, q_d9) {
  a <- p_d0 %>% transmute(slice_d0 = Slice_ID, mouse_d0 = mouse,
                          cells_d0 = cells, mol_d0 = mol, delta_d0 = delta)
  b <- q_d9 %>% transmute(slice_d9 = Slice_ID, mouse_d9 = mouse,
                          cells_d9 = cells, mol_d9 = mol, delta_d9 = delta)
  if (nrow(a) < MIN_SLICES || nrow(b) < MIN_SLICES) return(NULL)
  pr <- crossing(a, b) %>%
    mutate(shift = delta_d9 - delta_d0, same_mouse = mouse_d0 == mouse_d9)
  list(summary = tibble(
         n_slices_d0 = nrow(a), n_slices_d9 = nrow(b), n_pairs = nrow(pr),
         n_mice_d0 = n_distinct(a$mouse_d0), n_mice_d9 = n_distinct(b$mouse_d9),
         delta_d0 = median(a$delta_d0), delta_d9 = median(b$delta_d9),
         mean_shift = mean(pr$shift), median_shift = median(pr$shift),
         sd_shift = stats::sd(pr$shift), min_shift = min(pr$shift), max_shift = max(pr$shift),
         consistency = 100 * max(mean(pr$shift > 0), mean(pr$shift < 0))),
       pairs = pr)
}

# ---- one gene through one lineage, both designs -> 0-2 result rows --------------------
run_gene_lineage <- function(dat, gene, lg) {
  rt <- route_gene(dat, gene, lg)
  stub <- function(design, extra = list())
    tibble(gene = gene, lineage = lg$name, design = design, route = rt$route,
           marker = rt$marker, exclude_reason = rt$exclude_reason,
           p_elig = rt$p_elig %||% NA, p_avail = rt$p_avail %||% NA, p_med_mol = rt$p_med_mol %||% NA,
           d_elig = rt$d_elig %||% NA, d_avail = rt$d_avail %||% NA, d_med_mol = rt$d_med_mol %||% NA) %>%
      # NB: bind_cols(x, as_tibble(list())) returns ZERO rows - dplyr treats the empty
      # 0x0 tibble as a zero-row frame. Guard, or every excluded stub vanishes.
      { if (length(extra)) bind_cols(., as_tibble(extra)) else . }
  tag <- function(x, design, side = NA)
    if (!nrow(x)) x else x %>% mutate(gene = gene, lineage = lg$name, design = design,
                                      side = side, .before = 1)
  if (rt$route == "excluded")
    return(list(summary = stub("none"), pairs = tibble(), slices = tibble()))

  a0 <- delta_by_slice(dat, lg$parent,   gene) %>% filter(day == "Day 0")
  b9 <- delta_by_slice(dat, lg$daughter, gene) %>% filter(day == "Day 9")
  a9 <- delta_by_slice(dat, lg$parent,   gene) %>% filter(day == "Day 9")
  # slices dump keeps EVERY slice incl. ineligible ones (eligible flag says which)
  slc <- bind_rows(tag(a0, "any", "parent_D0"), tag(b9, "any", "daughter_D9"),
                   tag(a9, "any", "parent_D9"))
  p0 <- eligible_slices(a0); q9 <- eligible_slices(b9); p9 <- eligible_slices(a9)
  out <- list(); prs <- list()

  if (rt$route == "D1_full") {
    s1 <- pair_stats(p0, q9)
    if (!is.null(s1)) { out[[1]] <- stub("D1_lineage", s1$summary)
                        prs[[1]] <- tag(s1$pairs, "D1_lineage") }
    else out[[1]] <- stub("D1_lineage", list(status = "too_few_mice"))
    s2 <- pair_stats(p0, p9)
    if (!is.null(s2)) { out[[2]] <- stub("D2_persist", s2$summary)
                        prs[[2]] <- tag(s2$pairs, "D2_persist") }
  } else if (rt$route == "daughter_only") {
    # no Day 0 baseline exists (gene off in the healthy parent): report the within-
    # daughter position, to be judged against the blank null, not a shift.
    out[[1]] <- stub("daughter_only", list(
      n_slices_d9 = nrow(q9), n_mice_d9 = n_distinct(q9$mouse),
      delta_d9 = if (nrow(q9)) median(q9$delta) else NA_real_,
      consistency = if (nrow(q9)) 100 * max(mean(q9$delta > 0), mean(q9$delta < 0)) else NA_real_,
      status = "no_D0_baseline"))
  } else if (rt$route == "parent_only") {
    out[[1]] <- stub("parent_only", list(
      n_slices_d0 = nrow(p0), n_mice_d0 = n_distinct(p0$mouse),
      delta_d0 = if (nrow(p0)) median(p0$delta) else NA_real_,
      status = "lost_in_daughter"))
  }
  list(summary = bind_rows(out), pairs = bind_rows(prs), slices = slc)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- assemble the joined table for one layer ------------------------------------------
join_layer <- function(counts, meta, coords)
  counts %>% mutate(Cell_ID = as.character(Cell_ID)) %>%
    inner_join(meta,   by = c("Slice_ID", "Cell_ID")) %>%
    inner_join(coords, by = c("Slice_ID", "Cell_ID"))
