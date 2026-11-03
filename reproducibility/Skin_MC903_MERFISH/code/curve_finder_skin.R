# =====================================================================
# curve_finder_skin.R - the colon fitter, unchanged, with ONE override (the kidney twin).
# Sources morpho_improve/final/curve_finder.R verbatim (.connect, .project_all, .spike_score, .to_manual, the ladder).
# Override: fit_layer_xyt() selects the layer by EXACT label match (EPI / PAP / RET / SUB / SUBFAS; a prefix match
# would let "SUB" swallow "SUBFAS"). Everything else is the colon code.
# =====================================================================
.COLON_FITTER <- local({
  rel <- file.path("morpho_improve", "final", "curve_finder.R")
  roots <- c(getwd(), file.path(getwd(), ".."), file.path(getwd(), "..", ".."),
             local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
                     if (!length(a)) NULL else { d <- dirname(normalizePath(sub("^--file=", "", a[1]))); c(d, file.path(d, ".."), file.path(d, "..", "..")) } }))
  hit <- NULL
  for (r in roots) { p <- file.path(r, rel); if (!is.null(r) && file.exists(p)) { hit <- p; break } }
  if (is.null(hit)) stop("colon fitter (", rel, ") not found from: ", getwd())
  normalizePath(hit)
})
suppressMessages(source(.COLON_FITTER))

fit_layer_xyt <- function(slice_df, layer, min_cells = 30, force_loop = NULL, return_all = FALSE) {
  sel <- slice_df$Leiden_neigh == layer; sel[is.na(sel)] <- FALSE
  s <- slice_df[sel & is.finite(slice_df$x) & is.finite(slice_df$y), , drop = FALSE]
  keep <- names(which(table(s$Tier2) > min_cells)); s <- s[s$Tier2 %in% keep, , drop = FALSE]
  if (nrow(s) == 0) return(NULL)
  gcen <- c(stats::median(slice_df$x[is.finite(slice_df$x)]), stats::median(slice_df$y[is.finite(slice_df$y)]))
  cc  <- fit_candidates(as.matrix(s[, c("x","y")]), force_loop = force_loop, center = gcen)
  out <- .to_manual(cc$chosen, s, layer)
  if (return_all) { alts <- list(); for (nm in c("plain","reml","binned")) if (!is.null(cc[[nm]])) alts[[nm]] <- .to_manual(cc[[nm]], s, layer)
    attr(out, "candidates") <- alts; attr(out, "chosen_key") <- cc$chosen_key }
  out
}
