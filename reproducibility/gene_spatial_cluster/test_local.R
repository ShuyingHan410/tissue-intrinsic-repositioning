# LOCAL smoke test: one gene through the EXACT cluster code path (delta_core.R),
# using the already-extracted gene_counts_gt.csv instead of a cluster chunk.
#   Rscript test_local.R Wnt2b
# Output = the same schema results/part_NN.csv will contain, printed + saved.
args <- commandArgs(trailingOnly = TRUE)
GENE <- if (length(args)) args[1] else "Wnt2b"
Sys.setenv(OUT = ".")
source("delta_core.R")
COORD <- "../../output_threelayer/per_slice_morphogam"

gt <- read_csv("../gene_demo/gene_counts_gt.csv", show_col_types = FALSE, progress = FALSE) %>%
  mutate(Cell_ID = as.character(Cell_ID))
stopifnot(GENE %in% names(gt))
counts <- gt %>% select(Slice_ID, Cell_ID, all_of(GENE), blank_counts)
meta   <- gt %>% transmute(Slice_ID, Cell_ID, Tier1, Tier2,
                           day   = ifelse(grepl("_D9_", Slice_ID), "Day 9", "Day 0"),
                           mouse = sub("_slice_.*", "", Slice_ID),
                           total_counts)

cat(sprintf("== local smoke test: %s ==\n\n", GENE))
co_mu <- load_coords(COORD, "MU")
ot <- orientation_table(meta, co_mu)
cat(sprintf("orientation check (harmonised: lumen = 0, outer = 1): %d ok | %d INVERTED | %d not evaluable | %d slices\n\n",
            sum(ot$status == "ok"), sum(ot$status == "INVERTED"),
            sum(ot$status == "not_evaluable"), nrow(ot)))

joined <- list(MU = join_layer(counts, meta, co_mu), SM = {
  co_sm <- load_coords(COORD, "SM"); if (is.null(co_sm)) NULL else join_layer(counts, meta, co_sm) })

all <- lapply(LINEAGES, function(lg) {
  d <- joined[[lg$layer]]
  if (is.null(d)) { cat(sprintf("  [%s] no %s coordinates - skipped\n", lg$name, lg$layer))
                    return(list(summary=tibble(), pairs=tibble(), slices=tibble())) }
  run_gene_lineage(d, GENE, lg)
})
grab <- function(k) bind_rows(lapply(all, `[[`, k))
res <- grab("summary"); prs <- grab("pairs"); slc <- grab("slices")
cat("\n---- result rows (the schema every cluster part file uses) ----\n\n")
print(as.data.frame(res %>% mutate(across(where(is.numeric), ~round(., 4)))), row.names = FALSE)
write.csv(res, sprintf("results/local_test_%s.csv",        GENE), row.names = FALSE)
write.csv(prs, sprintf("results/local_test_%s_pairs.csv",  GENE), row.names = FALSE)
write.csv(slc, sprintf("results/local_test_%s_slices.csv", GENE), row.names = FALSE)
cat(sprintf("\n---- per-PAIR rows (%d) - use these to build any summary you like ----\n\n", nrow(prs)))
if (nrow(prs)) print(as.data.frame(prs %>% mutate(across(where(is.numeric), ~round(., 4))) %>% head(8)), row.names = FALSE)
cat(sprintf("\n---- per-SLICE deltas (%d) ----\n\n", nrow(slc)))
if (nrow(slc)) print(as.data.frame(slc %>% mutate(across(where(is.numeric), ~round(., 4))) %>% head(8)), row.names = FALSE)
types <- unique(unlist(lapply(LINEAGES, function(l) c(l$parent, l$daughter))))
cells <- bind_rows(lapply(names(joined), function(ly) {
  d <- joined[[ly]]; if (is.null(d)) return(tibble())
  d %>% filter(Tier2 %in% types) %>%
    transmute(layer = ly, Slice_ID, Cell_ID, Tier2, day, mouse, r,
              total_counts, count = .data[[GENE]]) }))
write.csv(cells, sprintf("results/local_test_%s_cells.csv", GENE), row.names = FALSE)
cat(sprintf("\n---- per-CELL rows (%d) - the rawest level ----\n\n", nrow(cells)))
print(as.data.frame(cells %>% filter(count > 0) %>% arrange(desc(count)) %>% head(5)), row.names = FALSE)
cat(sprintf("\nsaved -> results/local_test_%s{,_pairs,_slices,_cells}.csv\n", GENE))
