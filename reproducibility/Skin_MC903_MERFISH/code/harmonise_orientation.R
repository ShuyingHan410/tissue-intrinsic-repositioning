# =====================================================================
# harmonise_orientation.R - make r point the same way in every skin layer file (the kidney script, references swapped).
#   Rscript code/harmonise_orientation.R [output/morphogam_auto]
# CONVENTION: r = 0 on the SURFACE side, r = 1 on the DEEP side, in every compartment.
# Reference (geometry only): depth_ref stored in each file = distance to the keratinocyte sheet (PAP/RET/SUB/SUBFAS) or minus
# the distance to the nearest dermal cell (EPI). Spearman(r, depth_ref) must be positive; if not, r -> 1 - r. Report-only:
# in normal strips papillary fibroblasts should have smaller mean r than reticular within the dermis.
# =====================================================================
suppressMessages(library(data.table))
ROOT <- if (dir.exists("all_slices")) "." else ".."
FIT  <- if (length(commandArgs(TRUE))) file.path(ROOT, commandArgs(TRUE)[1]) else file.path(ROOT, "output", "morphogam_auto")
REVD <- sub("morphogam([^/]*)$", "morphogam\\1_review", FIT)
PAT  <- "^skin_per_slice_(.*)_(EPI|PAP|RET|SUB|SUBFAS)_xyt_mapped\\.csv$"
fs <- sort(list.files(FIT, pattern = "_xyt_mapped\\.csv$", full.names = TRUE)); cat(sprintf("%d layer files\n", length(fs)))
log <- list()
for (f in fs) {
  sid <- sub(PAT, "\\1", basename(f)); lay <- sub(PAT, "\\2", basename(f)); d <- fread(f)
  ref <- if ("depth_ref" %in% names(d)) d$depth_ref else d$depth
  rho <- suppressWarnings(cor(d$r_for_jsd, ref, method = "spearman")); flip <- isTRUE(rho < 0)
  if (flip) { d[, `:=`(r = 1 - r, r_original = 1 - r_original, r_for_jsd = 1 - r_for_jsd, morphogam_r_flip = TRUE)]; fwrite(d, f)
    for (rv in list.files(file.path(REVD, sid), pattern = sprintf("^%s_CHOSEN-.*_xyt_mapped\\.csv$", lay), full.names = TRUE)) fwrite(d, rv) }
  mr <- d[Tier2 %in% c("FIB_Papillary", "FIB_Reticular"), .(r = mean(r_for_jsd), n = .N), by = Tier2]
  pap <- mr[Tier2 == "FIB_Papillary" & n >= 30, r]; ret <- mr[Tier2 == "FIB_Reticular" & n >= 30, r]
  log[[length(log) + 1]] <- data.table(slice = sid, layer = lay, n = nrow(d), rho_before = round(rho, 3), rho_after = round(abs(rho), 3), flipped = flip,
                                       pap_mean_r = if (length(pap)) round(pap, 3) else NA_real_, ret_mean_r = if (length(ret)) round(ret, 3) else NA_real_,
                                       pap_lt_ret = if (length(pap) && length(ret)) pap < ret else NA)
  cat(sprintf("  %-20s %-6s %-8s rho %+.2f%s\n", sid, lay, if (flip) "FLIPPED" else "ok", abs(rho),
              if (length(pap) && length(ret)) sprintf("   papillary %.2f %s reticular %.2f", pap, if (pap < ret) "<" else ">", ret) else ""))
}
lg <- rbindlist(log); fwrite(lg, file.path(FIT, "_ORIENTATION_log.csv"))
cat(sprintf("\n%d of %d flipped; weak reference (|rho| < 0.3): %d; papillary < reticular where both present: %d of %d\n", sum(lg$flipped), nrow(lg),
            sum(lg$rho_after < 0.3), sum(lg$pap_lt_ret, na.rm = TRUE), sum(!is.na(lg$pap_lt_ret))))
