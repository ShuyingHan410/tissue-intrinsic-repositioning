# Harmonise r_for_jsd orientation across ALL per-slice coordinate files to the documented
# convention: inner (lumen) = 0 -> outer (serosa) = 1.
#
# Rule (geometric, works on every layer incl. SM/ME which have no epithelial anchors):
#   cor(r_for_jsd, distance from the layer centroid, spearman)
#     cor >  0.05  already inner=0/outer=1        -> leave
#     cor < -0.05  inverted (inner=1/outer=0)     -> r_for_jsd := 1 - r_for_jsd  (and r)
#     |cor| <= .05 ambiguous                      -> leave, flag in audit
# r_original is NEVER touched (provenance); morphogam_r_flip is set TRUE on flipped files.
# Backups: per_slice_morphogam*_preflip20260810. Audit -> orientation_harmonisation_audit.csv
suppressMessages({library(readr); library(dplyr)})
DIRS <- c("per_slice_morphogam", "per_slice_morphogam_NEW")
THR  <- 0.05
audit <- list()
for (D in DIRS) {
  fs <- list.files(D, pattern = "_(MU|SM|ME)_xyt_mapped\\.csv$", full.names = TRUE)
  for (f in fs) {
    d <- tryCatch(suppressMessages(read_csv(f, show_col_types = FALSE, progress = FALSE)),
                  error = function(e) NULL)
    row <- function(action, cor) tibble(dir = D, file = basename(f), cor_before = round(cor, 3),
                                        action = action)
    if (is.null(d) || !all(c("x","y","r_for_jsd") %in% names(d))) {
      audit[[length(audit)+1]] <- row("SKIP (unreadable/columns)", NA); next }
    ok <- is.finite(d$x) & is.finite(d$y) & is.finite(d$r_for_jsd)
    if (sum(ok) < 50) { audit[[length(audit)+1]] <- row("SKIP (<50 cells)", NA); next }
    cen  <- c(median(d$x[ok]), median(d$y[ok]))
    dist <- sqrt((d$x - cen[1])^2 + (d$y - cen[2])^2)
    cr <- suppressWarnings(cor(d$r_for_jsd[ok], dist[ok], method = "spearman"))
    if (!is.finite(cr)) { audit[[length(audit)+1]] <- row("SKIP (cor NA)", NA); next }
    if (cr < -THR) {
      d$r_for_jsd <- 1 - d$r_for_jsd
      if ("r" %in% names(d)) d$r <- 1 - d$r
      if ("morphogam_r_flip" %in% names(d)) d$morphogam_r_flip <- TRUE
      write_csv(d, f)
      audit[[length(audit)+1]] <- row("FLIPPED", cr)
    } else if (cr > THR) {
      audit[[length(audit)+1]] <- row("ok", cr)
    } else {
      audit[[length(audit)+1]] <- row("AMBIGUOUS (left as-is)", cr)
    }
  }
}
a <- bind_rows(audit)
write.csv(a, "orientation_harmonisation_audit.csv", row.names = FALSE)
cat("actions:\n"); print(table(a$dir, a$action))
cat("\nflipped files:\n")
print(a %>% filter(action == "FLIPPED") %>% arrange(dir, file) %>% as.data.frame(), row.names = FALSE)
cat("\nambiguous (review by eye):\n")
print(a %>% filter(grepl("AMBIG", action)) %>% as.data.frame(), row.names = FALSE)
