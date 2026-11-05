# Pool + score every method's Day0->Day9 partials three ways:
#   raw            comparator coupling MASS (as originally reported), ties -> min
#   colnorm        mass per unit of Day9 target mass (abundance-corrected), ties -> min
#   colnorm_avg    same, ties averaged (hard-matching methods produce ties)
# MorphoGAM (S) and the expression baseline are not couplings and are scored unchanged.
# Comparators ran on the adata grid (10 Day0 x 25 Day9 = 250 units); MorphoGAM/expression ran on
# the 21 Day9 slices with manual layers (210 units) -> comparators are restricted to those 21.
# Run from the benchmark root:  Rscript script/day0_day9/rescore_all_methods.R
suppressMessages({library(dplyr); library(readr); library(tibble); library(purrr); library(tidyr)})
source("script/common/recovery.R")
base <- "benchmark_out/day0_day9"
strip <- function(f) sub("[.]csv$", "", basename(f))
d9 <- sort(strip(list.files("all_slices", "_D9_.*[.]csv$", full.names = TRUE)))
d0 <- sort(strip(list.files("all_slices", "_D0_.*[.]csv$", full.names = TRUE)))
manual9 <- unique(sub("^DSS9_all_per_slice_", "", sub("_MU_xyt_mapped[.]csv$", "",
  list.files("output_threelayer/per_slice_morphogam", "^DSS9_all_.*_MU_xyt_mapped[.]csv$"))))
n9 <- length(d9); allu <- 0:(length(d0) * n9 - 1)
keep_units <- allu[sapply(allu, function(u) d9[u %% n9 + 1] %in% manual9)]
cat(sprintf("Day0 %d, Day9 adata %d, Day9 manual %d, comparator units kept %d of %d\n",
            length(d0), n9, length(manual9), length(keep_units), length(allu)))

read_parts <- function(m, keep = NULL) {
  fs <- sort(list.files(file.path(base, m, "partials"), pattern = "_Fibroblast_unit\\d+\\.csv$", full.names = TRUE))
  u <- as.integer(sub(".*unit(\\d+)\\.csv$", "\\1", fs)); if (!is.null(keep)) fs <- fs[u %in% keep]
  setNames(lapply(fs, function(f) as.matrix(read.csv(f, row.names = 1, check.names = FALSE))), basename(fs))
}
colnorm <- function(M) { cs <- colSums(M, na.rm = TRUE); cs[cs == 0] <- NA; sweep(M, 2, cs, "/") }
score_all <- function(mats, m, transform = identity, ties = "min") imap_dfr(mats, function(M, u) {
  M2 <- transform(M)
  pmap_dfr(GT_DAY0DAY9, function(origin, cognates, compartment) {
    if (!(origin %in% rownames(M2))) return(tibble())
    v <- M2[origin, ]; v <- v[is.finite(v)]; present <- intersect(unlist(cognates), names(v))
    if (length(v) < 2 || !length(present)) return(tibble())
    rk <- rank(-v, ties.method = ties); best <- present[which.min(rk[present])]
    tibble(method = m, unit = u, origin = origin, target = paste(unlist(cognates), collapse = "/"),
           rank = unname(rk[best]), n_cand = length(v), norm_rank = (rk[best] - 1) / (length(v) - 1))
  })
})
summ   <- function(s) s %>% group_by(method, origin, target) %>%
  summarise(n = n(), mean_norm = round(mean(norm_rank), 3), median_rank = median(rank),
            top1 = round(100 * mean(rank == 1)), top3 = round(100 * mean(rank <= 3)), .groups = "drop")
pooled <- function(s) s %>% group_by(method) %>%
  summarise(pairs = n(), mean_norm = round(mean(norm_rank), 3), top1 = round(100 * mean(rank == 1)),
            top3 = round(100 * mean(rank <= 3)), .groups = "drop")

# methods present on disk
own  <- intersect(c("morphogam", "expression"), list.dirs(base, recursive = FALSE, full.names = FALSE))
comp_all <- c("moscot", "moscot_ub", "paste2", "stalign", "slat", "cast_stack", "cast")
comp <- comp_all[sapply(comp_all, function(m) length(list.files(file.path(base, m, "partials"))) > 0)]
cat("own-scale methods:", own, "\ncomparators found:", comp, "\n")
own_m  <- lapply(own,  read_parts);             names(own_m)  <- own
comp_m <- lapply(comp, read_parts, keep_units); names(comp_m) <- comp
n_units <- sapply(comp_m, length); print(n_units)
if (any(n_units < length(keep_units))) warning("some comparator has fewer than ", length(keep_units), " kept units: incomplete run?")

s_raw    <- bind_rows(imap_dfr(own_m, score_all), imap_dfr(comp_m, function(mats, m) score_all(mats, m)))
s_col    <- bind_rows(imap_dfr(own_m, score_all), imap_dfr(comp_m, function(mats, m) score_all(mats, m, colnorm)))
s_colavg <- bind_rows(imap_dfr(own_m, function(mats, m) score_all(mats, m, identity, "average")),
                      imap_dfr(comp_m, function(mats, m) score_all(mats, m, colnorm, "average")))
options(width = 220)
cat("\n==== RAW: comparators on coupling mass (ties -> min)\n");            print(summ(s_raw) %>% arrange(origin, mean_norm)); print(pooled(s_raw) %>% arrange(mean_norm))
cat("\n==== ABUNDANCE-CORRECTED: column-normalised (ties -> min)\n");       print(summ(s_col) %>% arrange(origin, mean_norm)); print(pooled(s_col) %>% arrange(mean_norm))
cat("\n==== ABUNDANCE-CORRECTED, ties averaged\n");                        print(pooled(s_colavg) %>% arrange(mean_norm))
ties <- imap_dfr(comp_m, function(mats, m) tibble(method = m,
  pct_rows_tied_top_colnorm = round(100 * mean(sapply(mats, function(M) { M2 <- colnorm(M)
    v <- M2[rownames(M2) %in% c("Fibro 2", "Fibro 6", "Fibro 7"), , drop = FALSE]
    any(apply(v, 1, function(r) sum(r == max(r, na.rm = TRUE), na.rm = TRUE) > 1)) })))))
cat("\n==== share of kept pairs with a tie at the top of an origin row\n"); print(ties)
write_csv(bind_rows(s_raw %>% mutate(scoring = "raw"), s_col %>% mutate(scoring = "column-normalised"),
                    s_colavg %>% mutate(scoring = "column-normalised, ties averaged")),
          file.path(base, "per_pair_scores_all_methods.csv"))
write_csv(bind_rows(summ(s_raw) %>% mutate(scoring = "raw"), summ(s_col) %>% mutate(scoring = "column-normalised"),
                    summ(s_colavg) %>% mutate(scoring = "column-normalised, ties averaged")),
          file.path(base, "recovery_comparison_all_methods.csv"))
write_csv(bind_rows(pooled(s_raw) %>% mutate(scoring = "raw"), pooled(s_col) %>% mutate(scoring = "column-normalised"),
                    pooled(s_colavg) %>% mutate(scoring = "column-normalised, ties averaged")),
          file.path(base, "recovery_pooled_all_methods.csv"))
cat("\nwrote", file.path(base, c("per_pair_scores_all_methods.csv", "recovery_comparison_all_methods.csv", "recovery_pooled_all_methods.csv")), sep = "\n  ")
