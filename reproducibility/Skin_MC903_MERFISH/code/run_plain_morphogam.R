# =====================================================================
# run_plain_morphogam.R - PLAIN MorphoGAM on every skin strip x depth compartment (the kidney runner, constants swapped).
#
#   Rscript code/run_plain_morphogam.R                 all 16 strips, EPI PAP RET SUB SUBFAS
#   Rscript code/run_plain_morphogam.R s02_WT_Neck     one strip (regex on the slice id)
#   Rscript code/run_plain_morphogam.R s02 PAP         one strip, one layer
#   FIT_CAP=4000 Rscript code/run_plain_morphogam.R    change the subsample size
#
# "Plain" = ONE CurveFinder call per layer, no ladder, no REML, no binned fallback. Around it, only:
#   1. knn raised until the kNN graph is connected (.connect from the colon fitter);
#   2. layers above FIT_CAP cells fitted on a random subsample, every cell then projected (.project_all);
#   3. orientation: r = 0 at the skin SURFACE, r = 1 deep. The sign is set so that r increases with a landmark
#      depth reference (geometry, not a cell-type label): for the dermal / subcutaneous compartments the distance
#      to the keratinocyte sheet; for EPI itself that distance is ~0 for every cell, so EPI uses MINUS the distance
#      to the nearest dermal cell (deeper in the epidermis = closer to the dermis = larger r).
# Every strip is an open band -> loop = FALSE for all layers (no rings in skin).
#
# Output (own directory, like the kidney):
#   output/morphogam_auto/        skin_per_slice_<slice>_<LAYER>_xyt_mapped.csv  + _AUDIT_plain.csv
#   output/morphogam_auto_review/<slice>/<LAYER>_CHOSEN-plain_{curve,coordinate}_plot.png (+ CSV)
# Then:  Rscript -e 'rmarkdown::render("code/skin_gallery.Rmd")'
# =====================================================================
HERE <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd() })
setwd(normalizePath(file.path(HERE, "..")))                         # -> benchmark/Skin_MC903_MERFISH
suppressMessages(source(file.path("code", "curve_finder_skin.R")))  # MorphoGAM + the colon helpers
suppressMessages({ library(readr); library(ggplot2) })

args     <- commandArgs(TRUE)
SLICE_RX <- if (length(args) >= 1) args[1] else "."
LAYERS   <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("EPI", "PAP", "RET", "SUB", "SUBFAS")
FIT_CAP  <- as.integer(Sys.getenv("FIT_CAP", "5000"))
MIN_CELLS_LAYER <- 200                                              # a compartment thinner than this is not fitted
# optional env: SLICES (input dir) and OUT (output dir) so an alternative layer labelling can be run side by side, e.g.
#   SLICES=all_slices_B OUT=output/morphogam_B Rscript code/run_plain_morphogam.R s0 EPI,DER,SUB
IN_DIR <- Sys.getenv("SLICES", "all_slices"); OUT_DIR <- Sys.getenv("OUT", "output/morphogam_auto"); REV_DIR <- paste0(OUT_DIR, "_review")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE); dir.create(REV_DIR, showWarnings = FALSE, recursive = TRUE)
AUDIT <- file.path(OUT_DIR, "_AUDIT_plain.csv")
if (!file.exists(AUDIT)) writeLines("slice,time,layer,n,n_fit,knn,loop_arg,mode,spike,rho_depth,secs", AUDIT)

plain_fit <- function(xy, depth_ref, cap) {
  n <- nrow(xy); set.seed(1)
  fi  <- if (n > cap) sort(sample(n, cap)) else seq_len(n)
  xyf <- xy[fi, , drop = FALSE]
  k   <- .connect(xyf)$k
  fit <- CurveFinder(xyf, knn = k, loop = FALSE)                     # <- the one MorphoGAM call
  pr  <- .project_all(xy, as.numeric(fit$xyt$f1), as.numeric(fit$xyt$f2), as.numeric(fit$xyt$t), FALSE)
  rs  <- pr$rs; rho <- suppressWarnings(cor(rs, depth_ref, method = "spearman"))
  if (isTRUE(rho < 0)) { rs <- -rs; rho <- -rho }                   # surface = 0, deep = 1
  r01 <- (rs - min(rs)) / (max(rs) - min(rs) + 1e-12)
  out <- data.frame(x = xy[, 1], y = xy[, 2], t = pr$t, r = r01, r_original = r01, r_for_jsd = r01, ctr_x = pr$ctr_x, ctr_y = pr$ctr_y)
  attr(out, "mode") <- "open"; attr(out, "k_used") <- k; attr(out, "n_fit") <- nrow(xyf); attr(out, "rho") <- round(rho, 3)
  attr(out, "method") <- sprintf("plain-mg(open,knn=%d,nfit=%d)", k, nrow(xyf))
  attr(out, "spike")  <- round(.spike_score(xy, cbind(pr$ctr_x, pr$ctr_y)), 1)
  attr(out, "mg_plots") <- list(curve = fit$curve.plot, coordinate = fit$coordinate.plot)
  out
}

files <- grep(SLICE_RX, sort(list.files(IN_DIR, pattern = "^s\\d\\d_.*\\.csv$", full.names = TRUE)), value = TRUE)
if (!length(files)) stop("no strip matches '", SLICE_RX, "' in ", IN_DIR)
cat(sprintf("[plain MorphoGAM, skin] %d strip(s) x %s, cap %d cells\n", length(files), paste(LAYERS, collapse = "/"), FIT_CAP))

for (fl in files) {
  sid <- sub("\\.csv$", "", basename(fl))
  df  <- as.data.frame(suppressMessages(read_csv(fl, show_col_types = FALSE)))
  tp  <- as.character(df$Sample_type[1])
  sdir <- file.path(REV_DIR, sid); dir.create(sdir, showWarnings = FALSE, recursive = TRUE)
  for (lay in LAYERS) {
    out_csv <- file.path(OUT_DIR, sprintf("skin_per_slice_%s_%s_xyt_mapped.csv", sid, lay))
    if (file.exists(out_csv)) { cat("  skip (exists):", sid, lay, "\n"); next }
    s <- df[df$Leiden_neigh == lay & is.finite(df$x) & is.finite(df$y) & !is.na(df$Leiden_neigh), , drop = FALSE]
    s <- s[s$Tier2 %in% names(which(table(s$Tier2) > 30)), , drop = FALSE]     # the pipeline's >30-cell filter
    if (nrow(s) < MIN_CELLS_LAYER) { cat(sprintf("  %s %s: %d cells, not captured in this strip\n", sid, lay, nrow(s))); next }
    # depth reference for the orientation sign (see header): EPI -> minus distance to the nearest dermal cell
    if (lay == "EPI") { derm <- df[df$Leiden_neigh %in% c("PAP", "RET", "DER") & is.finite(df$x), c("x", "y")]
      s$depth_ref <- if (nrow(derm) >= 50) -FNN::get.knnx(as.matrix(derm), as.matrix(s[, c("x", "y")]), k = 1)$nn.dist[, 1] else s$depth
    } else s$depth_ref <- s$depth
    t0 <- Sys.time()
    f  <- tryCatch(plain_fit(as.matrix(s[, c("x", "y")]), s$depth_ref, FIT_CAP), error = function(e) { cat("  !!", sid, lay, "FAILED:", conditionMessage(e), "\n"); NULL })
    secs <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    if (is.null(f)) { write.table(data.frame(sid, tp, lay, nrow(s), NA, NA, FALSE, NA, NA, NA, secs), AUDIT, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE, quote = TRUE); next }
    out <- .to_manual(f, s, lay); out$morphogam_sample_label <- "skin"; out$depth <- s$depth; out$depth_ref <- s$depth_ref; out$layer_nn <- s$layer_nn
    write.csv(out, out_csv, row.names = FALSE)
    write.csv(out, file.path(sdir, sprintf("%s_CHOSEN-plain_xyt_mapped.csv", lay)), row.names = FALSE)
    mp <- attr(f, "mg_plots")
    ggsave(file.path(sdir, sprintf("%s_CHOSEN-plain_curve_plot.png", lay)), mp$curve, width = 7, height = 6, dpi = 90)
    set.seed(1); pidx <- if (nrow(out) > 8000) sort(sample(nrow(out), 8000)) else seq_len(nrow(out))
    gr <- ggplot(out[pidx, ][order(out$r_for_jsd[pidx]), ], aes(x, y, color = r_for_jsd)) + geom_point(size = .4) +
      scale_color_viridis_c(name = "r_for_jsd", limits = c(0, 1)) + coord_equal() + theme_bw() +
      ggtitle(sprintf("%s %s  plain MorphoGAM  r (surface=0 -> deep=1)", sid, lay))
    ggsave(file.path(sdir, sprintf("%s_CHOSEN-plain_coordinate_plot.png", lay)), gr, width = 7, height = 6, dpi = 90)
    write.table(data.frame(sid, tp, lay, nrow(s), attr(f, "n_fit"), attr(f, "k_used"), FALSE, attr(f, "mode"), attr(f, "spike"), attr(f, "rho"), secs),
                AUDIT, sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE, quote = TRUE)
    cat(sprintf("  %s %s: n=%d fit on %d, knn %d, spike %.1f, rho(r,depth) %.2f (%.0fs)\n", sid, lay, nrow(s), attr(f, "n_fit"), attr(f, "k_used"), attr(f, "spike"), attr(f, "rho"), secs))
  }
}
cat("\nDONE ->", OUT_DIR, "  audit:", AUDIT, "\n")
