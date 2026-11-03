# Shared Day0->Day9 fibroblast rank-recovery scorer.
# Applied IDENTICALLY to every method's per-pair matrices (rows = Day0 origin types,
# cols = Day9 candidate types) so MorphoGAM and the 4 alignment methods are compared
# on one footing. Ground truth = paper-supported polarizations; Fibro 6 counts as
# recovered if EITHER IAF 2 or IAF 3 ranks well (we take the better of the two).
suppressMessages({library(dplyr); library(readr); library(tibble); library(purrr)})

GT_DAY0DAY9 <- tibble::tribble(
  ~origin,   ~cognates,             ~compartment,
  "Fibro 2", list("IAF 4"),          "crypt top (MU)",
  "Fibro 6", list("IAF 2", "IAF 3"), "crypt base (MU)",   # best of the two (either polarization)
  "Fibro 7", list("Fibro 4"),        "submucosa (SM/ME)"
)

# row-wise rank of the BEST cognate among ALL finite candidates in row `origin`.
# rank 1 = best. Returns rank, the winning cognate, n candidates, normalized rank.
row_rank_best <- function(S, origin, cognates) {
  if (is.null(S) || !(origin %in% rownames(S))) return(tibble(rank=NA_real_, cognate=NA_character_, n_cand=NA_integer_))
  v <- S[origin, ]; v <- v[is.finite(v)]
  present <- intersect(unlist(cognates), names(v))
  if (!length(v) || !length(present)) return(tibble(rank=NA_real_, cognate=NA_character_, n_cand=length(v)))
  rk <- rank(-v, ties.method = "min")               # 1 = highest similarity
  best <- present[which.min(rk[present])]
  tibble(rank = unname(rk[best]), cognate = best, n_cand = length(v))
}

# score one per-pair matrix against the ground truth -> one row per GT linkage.
# `target` = the intended cognate label (identifies the linkage, so Fibro 6 -> IAF 2 and
# Fibro 6 -> IAF 3 are kept separate).
score_matrix <- function(S, gt = GT_DAY0DAY9) {
  purrr::pmap_dfr(gt, function(origin, cognates, compartment) {
    r <- row_rank_best(S, origin, cognates)
    tibble(origin = origin, target = paste(unlist(cognates), collapse = "/"),
           compartment = compartment,
           cognate = r$cognate, rank = r$rank, n_cand = r$n_cand,
           norm_rank = ifelse(is.finite(r$rank) & r$n_cand > 1,
                              (r$rank - 1) / (r$n_cand - 1), NA_real_))
  })
}

# read a method's per-pair partials (benchmark_out/day0_day9/<method>/partials/
# <method>_Fibroblast_unitNNN.csv), score each, tag with method + unit.
score_method_partials <- function(method, base_dir, lineage = "Fibroblast") {
  pdir <- file.path(base_dir, method, "partials")
  fs <- sort(list.files(pdir, pattern = sprintf("^%s_%s_unit\\d+\\.csv$", method, lineage),
                        full.names = TRUE))
  if (!length(fs)) return(tibble())
  purrr::imap_dfr(fs, function(f, i) {
    M <- tryCatch(as.matrix(read.csv(f, row.names = 1, check.names = FALSE)),
                  error = function(e) NULL)
    if (is.null(M)) return(tibble())
    score_matrix(M) %>% mutate(method = method, unit = basename(f), .before = 1)
  })
}

# pooled per-(method, origin, target) recovery summary across all testable pairs.
recovery_summary <- function(scored) {
  scored %>% filter(is.finite(rank)) %>%
    group_by(method, origin, target, compartment) %>%
    summarise(n_pairs = n(),
              cand_range = paste0(min(n_cand), "-", max(n_cand)),
              median_rank = median(rank),
              mean_norm_rank = round(mean(norm_rank, na.rm = TRUE), 3),
              pct_top1 = round(100 * mean(rank == 1)),
              pct_top3 = round(100 * mean(rank <= 3)),
              pct_top5 = round(100 * mean(rank <= 5)),
              .groups = "drop")
}
