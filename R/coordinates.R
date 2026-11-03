#' Tissue-intrinsic coordinates for one section and layer
#'
#' Fits a smooth medial curve through the cells of one anatomical layer with MorphoGAM
#' (`CurveFinder`) and returns, for every cell, the position along the curve `t` and the
#' normalized within-layer depth `r_for_jsd` in \eqn{[0, 1]} (0 = inner / surface side,
#' 1 = outer / deep side). The fitter wraps `CurveFinder` with a connectivity guard (the
#' smallest kNN that connects the graph), open/closed detection (`loop = "auto"`), REML
#' re-smoothing of the medial curve with projection of all cells, and a binned-centreline
#' fallback for layers the smooth fit cannot trace. This is the automatic fitter; the
#' colon coordinates in the paper were produced with `CurveFinder` directly, escalating to
#' a larger `knn` and to `CurveFinderInteractive` where the default fit failed
#' (see `reproducibility/three_layer.Rmd`).
#'
#' @param slice_df data.frame with one row per cell of the whole section and at least the
#'   columns `x`, `y`, the layer label column and a population column (`Tier2`);
#'   `Slice_ID`, `Cell_ID`, `Tier1`, `Tier3`, `Sample_type`, `Mouse_ID` are carried over
#'   when present.
#' @param layer the layer to fit.
#' @param layer_col name of the column holding the layer label.
#' @param match `"prefix"` matches labels that start with `layer` (colon neighbourhood
#'   labels such as `MU_1`), `"exact"` requires equality (skin compartments).
#' @param min_cells populations with at most this many cells in the layer are excluded
#'   before fitting.
#' @param force_loop `NULL` lets MorphoGAM decide between an open and a closed curve;
#'   `TRUE`/`FALSE` forces it.
#' @return a data.frame in the per-section coordinate format used throughout the
#'   repository (`x, y, t, r, f1, f2, ..., r_original, r_for_jsd, morphogam_*`), with the
#'   fit provenance in its attributes (`method`, `mode`, `quality`, `spike`, `flag`).
#'   `NULL` when the layer has no population passing `min_cells`.
#' @examples
#' \dontrun{
#' cells <- read.csv("all_slices/062921_D0_m3a_1_slice_1.csv")
#' xyt <- fit_layer_coordinates(cells, "MU", layer_col = "Leiden_neigh")
#' head(xyt[, c("Cell_ID", "Tier2", "t", "r_for_jsd")])
#' }
#' @export
fit_layer_coordinates <- function(slice_df, layer, layer_col = "Leiden_neigh",
                                  match = c("prefix", "exact"), min_cells = 30,
                                  force_loop = NULL) {
  match <- match.arg(match)
  stopifnot(all(c("x", "y", layer_col, "Tier2") %in% names(slice_df)))
  lab <- as.character(slice_df[[layer_col]])
  sel <- if (match == "prefix") startsWith(lab, layer) else lab == layer
  sel[is.na(sel)] <- FALSE
  s <- slice_df[sel & is.finite(slice_df$x) & is.finite(slice_df$y), , drop = FALSE]
  keep <- names(which(table(s$Tier2) > min_cells))
  s <- s[s$Tier2 %in% keep, , drop = FALSE]
  if (nrow(s) == 0) return(NULL)
  gx <- slice_df$x[is.finite(slice_df$x)]; gy <- slice_df$y[is.finite(slice_df$y)]
  gcen <- c(stats::median(gx), stats::median(gy))
  cc <- fit_candidates(as.matrix(s[, c("x", "y")]), force_loop = force_loop, center = gcen)
  .to_manual(cc$chosen, s, layer)
}

#' Harmonise the orientation of a within-layer coordinate
#'
#' Flips `r` to `1 - r` when it decreases with distance from the section centre, so that
#' 0 is the inner (lumen / surface) side and 1 the outer (deep) side in every section.
#' This is the rule applied to all colon coordinate files
#' (`reproducibility/output_threelayer/flip_orientation.R`).
#'
#' @param r numeric vector of positions in \eqn{[0, 1]}.
#' @param x,y cell coordinates.
#' @param centre optional centre of the whole section (defaults to the median of `x` and
#'   `y` of the cells supplied).
#' @param threshold Spearman correlations between `-threshold` and `threshold` are
#'   treated as ambiguous and left unchanged.
#' @return a list with `r` (possibly flipped), `flipped` (logical) and `rho`.
#' @export
orient_outward <- function(r, x, y, centre = NULL, threshold = 0.05) {
  if (is.null(centre)) centre <- c(stats::median(x, na.rm = TRUE), stats::median(y, na.rm = TRUE))
  d <- sqrt((x - centre[1])^2 + (y - centre[2])^2)
  rho <- suppressWarnings(stats::cor(r, d, method = "spearman", use = "complete.obs"))
  flip <- is.finite(rho) && rho < -threshold
  list(r = if (flip) 1 - r else r, flipped = flip, rho = rho)
}
