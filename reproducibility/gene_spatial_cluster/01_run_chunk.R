# Step 1 (array task): all genes in one chunk x all lineages x both designs.
#   Rscript 01_run_chunk.R <chunk_id>        -> results/part_<id>.csv
args <- commandArgs(trailingOnly = TRUE)
cid  <- as.integer(args[1]); stopifnot(is.finite(cid))
OUT  <- Sys.getenv("OUT", "."); COORD <- Sys.getenv("COORD_DIR", "../../output_threelayer/per_slice_morphogam")
source(file.path(OUT, "delta_core.R"))

meta   <- read_csv(file.path(OUT, "meta.csv.gz"), show_col_types = FALSE, progress = FALSE) %>%
  mutate(Cell_ID = as.character(Cell_ID))
cf <- file.path(OUT, "chunks", sprintf("chunk_%02d.csv.gz", cid))
if (!file.exists(cf)) stop(sprintf(
  "chunk_%02d.csv.gz does not exist - array size exceeds the %d chunks written by 00_extract_all.py",
  cid, length(list.files(file.path(OUT, "chunks"), "^chunk_.*csv.gz$"))))
counts <- read_csv(cf, show_col_types = FALSE, progress = FALSE) %>%
  mutate(Cell_ID = as.character(Cell_ID))
genes  <- setdiff(names(counts), c("Slice_ID", "Cell_ID"))
cat(sprintf("chunk %d: %d genes\n", cid, length(genes)))

layers <- unique(vapply(LINEAGES, function(l) l$layer, ""))
joined <- setNames(lapply(layers, function(ly) {
  co <- load_coords(COORD, ly)
  if (is.null(co)) return(NULL)
  join_layer(counts, meta, co)
}), layers)

# orientation audit (MU only; written once per chunk, identical across chunks)
if (cid == 1 && !is.null(joined[["MU"]])) {
  ot <- orientation_table(meta, load_coords(COORD, "MU"))
  write.csv(ot, file.path(OUT, "results", "orientation_check.csv"), row.names = FALSE)
  cat(sprintf("orientation: %d ok | %d INVERTED | %d not evaluable (no epithelium in MU) | %d slices\n",
              sum(ot$status == "ok"), sum(ot$status == "INVERTED"),
              sum(ot$status == "not_evaluable"), nrow(ot)))
  if (any(ot$status == "INVERTED")) stop(
    "INVERTED slices detected - every delta in them would carry the wrong sign:\n  ",
    paste(ot$Slice_ID[ot$status == "INVERTED"], collapse = "\n  "))
}

all <- lapply(genes, function(gn) lapply(LINEAGES, function(lg) {
  d <- joined[[lg$layer]]
  if (is.null(d)) return(list(summary = tibble(), pairs = tibble(), slices = tibble()))
  run_gene_lineage(d, gn, lg)
}))
grab <- function(k) bind_rows(lapply(all, function(g) bind_rows(lapply(g, `[[`, k))))
res <- grab("summary"); prs <- grab("pairs"); slc <- grab("slices")
dir.create(file.path(OUT, "results"), showWarnings = FALSE)
write.csv(res, file.path(OUT, "results", sprintf("part_%02d.csv",   cid)), row.names = FALSE)
write.csv(prs, file.path(OUT, "results", sprintf("pairs_%02d.csv",  cid)), row.names = FALSE)
write.csv(slc, file.path(OUT, "results", sprintf("slices_%02d.csv", cid)), row.names = FALSE)
cat(sprintf("chunk %02d: %d summary rows | %d pair rows | %d slice rows\n",
            cid, nrow(res), nrow(prs), nrow(slc)))

# ---- OPTIONAL per-CELL dump: DUMP_CELLS="Wnt2b,Mmp10" (or "all") ---------------------
dc <- Sys.getenv("DUMP_CELLS", "")
if (nzchar(dc)) {
  want <- if (dc == "all") genes else intersect(trimws(strsplit(dc, ",")[[1]]), genes)
  types <- unique(unlist(lapply(LINEAGES, function(l) c(l$parent, l$daughter))))
  cells <- bind_rows(lapply(names(joined), function(ly) {
    d <- joined[[ly]]; if (is.null(d)) return(tibble())
    d %>% filter(Tier2 %in% types) %>%
      select(Slice_ID, Cell_ID, Tier2, day, mouse, r, total_counts, all_of(want)) %>%
      mutate(layer = ly, .before = 1)
  }))
  if (nrow(cells)) {
    write.csv(cells, file.path(OUT, "results", sprintf("cells_%02d.csv", cid)), row.names = FALSE)
    cat(sprintf("           + %d per-CELL rows for %s -> results/cells_%02d.csv\n",
                nrow(cells), paste(want, collapse=","), cid))
  }
}
