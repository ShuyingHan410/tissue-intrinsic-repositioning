# =====================================================================
# export_slices.R - data/mc903_cells.csv.gz  ->  Skin_MC903_MERFISH/all_slices/<strip>.csv
# =====================================================================
# Same per-slice column schema the colon pipeline reads (the skin twin of Kidney_Xenium/code/export_slices.R),
# so morpho_improve/final/curve_finder.R and the plain-MorphoGAM runner work on this tissue unchanged.
#
# Column mapping (colon name  <-  skin field):
#   Leiden_neigh  <-  layer_depth (EPI / PAP / RET / SUB / SUBFAS)   the LAYER the fitter works on; built from
#                                                                     depth landmarks only (code/skin_layers.R)
#   Tier2         <-  celltype (42 fine types)                        what the >30-cell filter counts
#   Tier1         <-  broad class (Keratinocyte / Fibroblast / Immune / ...)
#   Tier3         <-  celltype_broad (26)
#   Slice_ID      <-  s<rank>_<sample>  rank-first so alphabetical order is normal -> MC903 time course -> OXA
#   Cell_ID       <-  cell id;  Mouse_ID <- sample (one section per strip, mouse unknown)
#   Sample_type   <-  timepoint;  Technical_repeat_number <- 1
#   extra: depth, layer_nn, cc11, n_counts, site, condition, rundate
#
#   Rscript code/export_slices.R        (run from benchmark/Skin_MC903_MERFISH)
# =====================================================================
suppressMessages(library(data.table))
HERE <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
setwd(normalizePath(file.path(HERE, "..")))            # -> benchmark/Skin_MC903_MERFISH
source("code/skin_layers.R")

IN <- "data/mc903_cells.csv.gz"; OUT <- "all_slices"
if (!file.exists(IN)) stop("not found: ", IN, " - run code/extract_cells.R first")
dir.create(OUT, showWarnings = FALSE)

ORDER <- c("WT_Ear","WT_Neck","WT_Flank","d0_WT_neck","d3_MC903_neck","d6_MC903_neck","d9_MC903_Ear","d9_MC903_Neck","d9_MC903_Flank",
           "d12_MC903_neck","d10_OXA_Ear","d10_OXA_Neck","d10_OXA_Flank","d16_OXA_Ear","d16_OXA_Neck","d16_OXA_Flank")
d <- fread(IN)
d[, class := fcase(grepl("^KC", celltype_broad), "Keratinocyte", grepl("^HF|^SG", celltype_broad), "HairFollicle_Sebaceous",
                   grepl("^FIB", celltype_broad), "Fibroblast", celltype_broad == "Adipocyte", "Adipocyte", celltype_broad == "Chondrocyte", "Cartilage",
                   celltype_broad %in% c("SM","Pericyte","VEC","LEC","Schwann","MEL"), "Vessel_Muscle_Nerve", default = "Immune")]
d <- assign_skin_layers(d)
B <- attr(d, "bounds"); fwrite(B[order(match(sample, ORDER))], file.path(OUT, "_layer_bounds.csv"))
d[, Slice_ID := sprintf("s%02d_%s", match(sample, ORDER), sample)]

cat(sprintf("%s cells, %d strips, layers: %s\n", format(nrow(d), big.mark = ","), uniqueN(d$Slice_ID), paste(LAY5, collapse = ", ")))
for (sid in sort(unique(d$Slice_ID))) {
  s <- d[Slice_ID == sid, .(x, y, Leiden_neigh = as.character(layer_depth), Tier1 = class, Tier2 = celltype, Tier3 = celltype_broad,
                            Slice_ID, Cell_ID = cell, Mouse_ID = sample, Sample_type = timepoint, Technical_repeat_number = 1L,
                            depth = round(depth, 1), layer_nn = as.character(layer_nn), cc11, n_counts, site, condition, rundate)]
  fwrite(s, file.path(OUT, paste0(sid, ".csv")))
  tb <- table(s$Leiden_neigh)[LAY5]; tb[is.na(tb)] <- 0
  cat(sprintf("  %-22s %6d cells | %s\n", sid, nrow(s), paste(sprintf("%s=%d", LAY5, tb), collapse = " ")))
}
cat(sprintf("\nwrote %d files to %s/ (+ _layer_bounds.csv)\n", uniqueN(d$Slice_ID), normalizePath(OUT)))
