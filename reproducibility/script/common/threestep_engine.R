# Shared three-step spatial-similarity engine, used by BOTH benchmark designs:
#   designA_halfsplit/threestep.Rmd  (one slice -> stratified halves A/B)
#   designB_slicepair/threestep.Rmd  (two full slices -> A = slice i, B = slice j)
# Only the DATA PREP differs (make_split vs make_pair) + the run loop + the type
# orders; the engine below (per-layer KDE/JSD -> weighted integration -> layer-
# overlap correction) is identical, so both designs score your method the same way.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tibble)
  library(tidyr); library(purrr); library(ggplot2); library(knitr); library(pheatmap)
})

# Robust path to the per-slice MorphoGAM outputs (works from any working dir /
# on the cluster). Add your cluster path here if it differs.
per_slice_xyt_dir <- c(
  Sys.getenv("COORD_DIR", ""),               # explicit override
  "output_threelayer/per_slice_morphogam",
  "../output_threelayer/per_slice_morphogam",
  "../../output_threelayer/per_slice_morphogam"
) %>% { .[which(nzchar(.) & dir.exists(.))[1]] }
stopifnot(!is.na(per_slice_xyt_dir))

layers <- c("MU", "SM", "ME")

# Tuned parameters (final method) -- shared by both designs.
CFG <- list(
  min_cells      = 30,   # >=30 cells per type per layer in BOTH sides
  bw_mult        = 2,    # KDE bandwidth = bw.nrd0(pooled) * 2
  bw_floor       = 0.05,
  grid_n         = 512,
  sigma_quantile = 0.5,  # S = exp(-D/sigma), sigma = median(D); rank-neutral
  gamma          = 1     # layer-overlap exponent
)

# read one slice's cells for a lineage (Tier1 matched by regex)
read_lineage <- function(slice_id, sample_prefix, lineage_regex) {
  bind_rows(lapply(layers, function(ly) {
    f <- file.path(per_slice_xyt_dir,
                   paste0(sample_prefix, "_per_slice_", slice_id, "_", ly,
                          "_xyt_mapped.csv"))
    if (!file.exists(f)) return(tibble())
    read_csv(f, show_col_types = FALSE)
  })) %>%
    mutate(Slice_ID = as.character(Slice_ID), layer = as.character(layer),
           Tier1 = as.character(Tier1), Tier2 = as.character(Tier2),
           r_for_jsd = as.numeric(r_for_jsd)) %>%
    filter(Slice_ID == slice_id, layer %in% layers,
           grepl(lineage_regex, Tier1, ignore.case = TRUE),
           is.finite(r_for_jsd), r_for_jsd >= 0, r_for_jsd <= 1)
}

# Jensen-Shannon distance between two samples via Gaussian KDE
js_dist_kde <- function(x, y, bw, grid_n = 512) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  dx <- density(x, bw = bw, from = 0, to = 1, n = grid_n)$y
  dy <- density(y, bw = bw, from = 0, to = 1, n = grid_n)$y
  px <- dx / sum(dx); py <- dy / sum(dy); m <- 0.5 * (px + py)
  kl <- function(p, q) { i <- p > 0 & q > 0; sum(p[i] * log(p[i] / q[i])) }
  sqrt(0.5 * kl(px, m) + 0.5 * kl(py, m))
}

order_by <- function(mat, ord) {
  ro <- c(intersect(ord, rownames(mat)), setdiff(rownames(mat), ord))
  co <- c(intersect(ord, colnames(mat)), setdiff(colnames(mat), ord))
  mat[ro, co, drop = FALSE]
}

# THE THREE-STEP ENGINE. `dat` must have columns: layer, Tier2, value, grp (A/B).
three_step <- function(dat, cfg, type_order) {
  per_layer <- function(ln) {
    dA <- dat %>% filter(layer == ln, grp == "A")
    dB <- dat %>% filter(layer == ln, grp == "B")
    kA <- dA %>% count(Tier2) %>% filter(n >= cfg$min_cells)
    kB <- dB %>% count(Tier2) %>% filter(n >= cfg$min_cells)
    dA <- dA %>% semi_join(kA, by = "Tier2"); dB <- dB %>% semi_join(kB, by = "Tier2")
    rs <- sort(unique(dA$Tier2)); cs <- sort(unique(dB$Tier2))
    D <- matrix(NA_real_, length(rs), length(cs), dimnames = list(rs, cs))
    if (!length(rs) || !length(cs)) return(list(S = D))
    rA <- split(dA$value, dA$Tier2); rB <- split(dB$value, dB$Tier2)
    bw <- bw.nrd0(c(dA$value, dB$value)) * cfg$bw_mult
    if (!is.finite(bw) || bw <= 0) bw <- cfg$bw_floor
    for (i in rs) for (j in cs) D[i, j] <- js_dist_kde(rA[[i]], rB[[j]], bw, cfg$grid_n)
    sig <- as.numeric(quantile(D, cfg$sigma_quantile, na.rm = TRUE))
    if (!is.finite(sig) || sig <= 0) sig <- 1
    list(S = exp(-D / sig))
  }
  J <- setNames(lapply(layers, per_layer), layers)
  ar <- sort(unique(unlist(lapply(J, function(x) rownames(x$S)))))
  ac <- sort(unique(unlist(lapply(J, function(x) colnames(x$S)))))
  if (!length(ar) || !length(ac)) return(NULL)

  Ac <- dat %>% filter(grp == "A") %>% count(layer, Tier2, name = "n")
  Bc <- dat %>% filter(grp == "B") %>% count(layer, Tier2, name = "n")
  z <- function() matrix(0, length(ar), length(ac), dimnames = list(ar, ac))
  Sn <- z(); Sd <- z()
  for (ly in layers) {
    S <- J[[ly]]$S; nA <- Ac %>% filter(layer == ly); nB <- Bc %>% filter(layer == ly)
    if (!length(rownames(S)) || !length(colnames(S))) next
    for (i in rownames(S)) for (j in colnames(S)) {
      s <- S[i, j]; if (!is.finite(s)) next
      ni <- nA %>% filter(Tier2 == i) %>% pull(n); nj <- nB %>% filter(Tier2 == j) %>% pull(n)
      if (!length(ni) || !length(nj)) next
      w <- sqrt(ni * nj); if (!is.finite(w) || w <= 0) next
      Sn[i, j] <- Sn[i, j] + w * s; Sd[i, j] <- Sd[i, j] + w
    }
  }
  Ss <- Sn / Sd; Ss[Sd == 0] <- NA

  lp <- function(cd) cd %>% group_by(Tier2) %>% mutate(p = n / sum(n)) %>% ungroup()
  Alp <- lp(Ac); Blp <- lp(Bc)
  LO <- matrix(NA_real_, nrow(Ss), ncol(Ss), dimnames = dimnames(Ss))
  for (i in rownames(Ss)) for (j in colnames(Ss)) {
    pi <- Alp %>% filter(Tier2 == i) %>% select(layer, pi = p)
    pj <- Blp %>% filter(Tier2 == j) %>% select(layer, pj = p)
    pp <- full_join(pi, pj, by = "layer") %>% mutate(pi = replace_na(pi, 0), pj = replace_na(pj, 0))
    LO[i, j] <- sum(sqrt(pp$pi * pp$pj))
  }
  order_by(Ss * (LO ^ cfg$gamma), type_order)
}

row_rank <- function(S, q, cand) {
  if (is.null(S) || !(q %in% rownames(S))) return(c(rank = NA, sim = NA, n = NA))
  v <- S[q, ]; v <- v[is.finite(v)]
  if (!(cand %in% names(v))) return(c(rank = NA, sim = NA, n = length(v)))
  c(rank = unname(rank(-v, ties.method = "min")[cand]), sim = unname(v[cand]), n = length(v))
}

# Average a list of (possibly differently sized) matrices over a fixed type order.
mean_matrix <- function(mats, type_order) {
  rows <- intersect(type_order, unique(unlist(lapply(mats, rownames))))
  cols <- intersect(type_order, unique(unlist(lapply(mats, colnames))))
  acc <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols)); cnt <- acc
  for (S in mats) for (i in intersect(rows, rownames(S))) for (j in intersect(cols, colnames(S))) {
    v <- S[i, j]; if (is.finite(v)) { acc[i, j] <- acc[i, j] + v; cnt[i, j] <- cnt[i, j] + 1 }
  }
  M <- acc / cnt; M[cnt == 0] <- NA; M
}

# Summary table for the expected pairs across units (slices or slice pairs).
pair_rank_table <- function(mats, pairs) {
  pmap_dfr(pairs, function(query, candidate, relationship) {
    rk <- map_dfr(mats, function(S) {
      r <- row_rank(S, query, candidate)
      tibble(rank = unname(r["rank"]), sim = unname(r["sim"]), n = unname(r["n"]))
    }) %>% filter(is.finite(rank))
    rk <- rk %>% mutate(norm_rank = ifelse(n > 1, (rank - 1) / (n - 1), 0))
    tibble(
      pair = paste0(query, " -> ", candidate), relationship = relationship,
      n_units = nrow(rk), cand_range = paste0(min(rk$n), "-", max(rk$n)),
      mean_rank = round(mean(rk$rank), 2), norm_rank = round(mean(rk$norm_rank), 3),
      pct_top1 = round(100 * mean(rk$rank == 1)), pct_top3 = round(100 * mean(rk$rank <= 3)),
      mean_similarity = round(mean(rk$sim), 3),
      rank_band = case_when(mean(rk$norm_rank) <= 0.10 ~ "rank 1 (identity)",
                            mean(rk$norm_rank) <= 0.50 ~ "top half", TRUE ~ "bottom half"))
  })
}

# self-match rank-1 rate over the diagonal (comparable to the OT methods' metric)
diag_rank1 <- function(M) {
  ok <- 0; tot <- 0
  for (i in rownames(M)) if (i %in% colnames(M)) {
    v <- M[i, ]; v <- v[is.finite(v)]
    if (i %in% names(v)) { tot <- tot + 1; if (rank(-v, ties.method = "min")[i] == 1) ok <- ok + 1 }
  }
  if (tot == 0) NA_real_ else ok / tot
}

plot_lineage_heatmap <- function(M, title) {
  rank_lab <- t(apply(M, 1, function(v) {
    out <- rep("", length(v)); fin <- is.finite(v)
    out[fin] <- as.character(rank(-v[fin], ties.method = "min")); out
  }))
  dimnames(rank_lab) <- dimnames(M)
  pheatmap(M, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("white", "gold", "red"))(100),
           display_numbers = rank_lab, fontsize_number = 9, na_col = "grey92",
           fontsize_row = 9, fontsize_col = 9, angle_col = 45, border_color = "white",
           main = title)
}
