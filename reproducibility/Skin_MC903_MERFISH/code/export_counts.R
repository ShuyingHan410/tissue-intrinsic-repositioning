# Export the 500-gene MERFISH counts of every MERFISH cell (joins to all_slices/<strip>.csv on Cell_ID).
#   Rscript code/export_counts.R        (run from Skin_MC903_MERFISH/)
suppressMessages({ library(SeuratObject); library(Matrix); library(data.table) })
o <- readRDS("data/MERFISH_sc_int_Normal_MC903_OXA_annotated.rds")
mi <- which(o$assay == "MERFISH"); cnt <- o@assays$RNA@layers$counts[, mi]       # v5 layers carry no dimnames
idx <- which(rowSums(cnt > 0) > 0); panel <- rownames(o@assays$RNA)[idx]; cnt <- cnt[idx, ]
cat("MERFISH cells:", length(mi), " panel genes:", length(panel), "\n")
m <- as.data.table(as.matrix(t(cnt))); setnames(m, panel); m[, Cell_ID := colnames(o)[mi]]; setcolorder(m, "Cell_ID")
m[, total_counts := rowSums(.SD), .SDcols = panel]
fwrite(m, "data/mc903_counts.csv.gz"); writeLines(panel, "data/panel_genes.txt"); cat("wrote data/mc903_counts.csv.gz\n")
