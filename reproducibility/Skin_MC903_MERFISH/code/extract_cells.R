# Export per-cell metadata + spatial coordinates from the Imanishi 2025 Seurat object (Mendeley s3hgv4fp7c).
#   Rscript code/extract_cells.R            (run from Skin_MC903_MERFISH/)
suppressMessages({ library(SeuratObject); library(data.table) })
ROOT <- if (dir.exists("data")) "." else ".."
obj <- readRDS(file.path(ROOT, "data", "MERFISH_sc_int_Normal_MC903_OXA_annotated.rds"))
cat("class:", class(obj), " cells:", ncol(obj), " features:", nrow(obj), "\n")
cat("assays:", paste(names(obj@assays), collapse = ", "), "\n")
cat("images:", paste(names(obj@images), collapse = ", "), "\n")
md <- as.data.table(obj@meta.data, keep.rownames = "cell")
cat("\nmeta columns:\n"); print(names(md))
for (cn in names(md)) {
  v <- md[[cn]]
  if (is.character(v) || is.factor(v) || is.logical(v)) {
    u <- length(unique(v)); if (u <= 60) { cat(sprintf("\n[%s] %d levels:\n", cn, u)); print(sort(table(v, useNA = "ifany"), decreasing = TRUE)) }
    else cat(sprintf("\n[%s] %d levels, e.g. %s\n", cn, u, paste(head(unique(v), 4), collapse = " | ")))
  }
}
# coordinates live in the metadata (center_x / center_y, per-run slide frame); no FOV images in the object
cat("\n[cellcharter_cluster.k_11]\n"); print(table(md$cellcharter_cluster.k_11, useNA = "ifany"))
cat("\nsample x timepoint x condition x rundate (MERFISH only):\n")
print(md[assay == "MERFISH", .N, by = .(sample, timepoint, condition, rundate)][order(rundate, timepoint, sample)])
out <- md[assay == "MERFISH"]; out[, `:=`(center_x = as.numeric(center_x), center_y = as.numeric(center_y), volume = as.numeric(volume))]
cat("\nextent per sample (um):\n"); print(out[, .(n = .N, xmin = round(min(center_x)), xmax = round(max(center_x)), ymin = round(min(center_y)), ymax = round(max(center_y))), by = sample][order(sample)])
fwrite(out[, .(cell, sample, timepoint, condition, rundate, site = sub(".*_", "", sample), celltype_broad = intres1.1celltypebroad,
              celltype = celltype_subtype, cc11 = cellcharter_cluster.k_11, n_counts = nCount_RNA, n_genes = nFeature_RNA, volume,
              x = center_x, y = center_y)], file.path(ROOT, "data", "mc903_cells.csv.gz"))
cat("wrote", file.path(ROOT, "data", "mc903_cells.csv.gz"), nrow(out), "MERFISH cells\n")
