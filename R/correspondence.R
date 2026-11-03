#' Jensen–Shannon distance between two positional samples
#'
#' Gaussian kernel density estimates of two samples on \eqn{[0, 1]} are compared by the
#' Jensen–Shannon distance (square root of the Jensen–Shannon divergence, natural log).
#'
#' @param x,y numeric vectors of positions in \eqn{[0, 1]}.
#' @param bw kernel bandwidth.
#' @param grid_n number of grid points.
#' @return a non-negative scalar, `NA` if either sample has fewer than two values.
#' @export
js_distance <- function(x, y, bw, grid_n = 512) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  dx <- stats::density(x, bw = bw, from = 0, to = 1, n = grid_n)$y
  dy <- stats::density(y, bw = bw, from = 0, to = 1, n = grid_n)$y
  px <- dx / sum(dx); py <- dy / sum(dy); m <- 0.5 * (px + py)
  kl <- function(p, q) { i <- p > 0 & q > 0; sum(p[i] * log(p[i] / q[i])) }
  sqrt(0.5 * kl(px, m) + 0.5 * kl(py, m))
}

#' Anatomical similarity between cell populations of two sections or conditions
#'
#' Implements the population-correspondence score of the paper. Within each anatomical
#' layer, positional distributions of every population of group A and every population of
#' group B are compared by [js_distance()] (bandwidth `bw_mult` times the Silverman
#' bandwidth of the pooled positions); distances become similarities
#' \eqn{S_h = \exp(-D_h/\sigma_h)} with \eqn{\sigma_h} the `sigma_quantile` quantile of
#' the layer's distances. Layers are combined with weights \eqn{w_h = \sqrt{N_a^h N_b^h}}
#' and the result is multiplied by the layer-occupancy overlap
#' \eqn{O_{ab} = \sum_h \sqrt{\pi_a^h \pi_b^h}} raised to `gamma`.
#'
#' @param data data.frame with one row per cell and columns for the layer, the
#'   population label, the position and the group.
#' @param groups character vector of length 2 naming the two groups to compare (rows of
#'   the result are populations of `groups[1]`, columns those of `groups[2]`). Defaults to
#'   the first two distinct values of the group column.
#' @param layer,population,r,group column names.
#' @param min_cells populations with fewer cells than this in a layer are dropped from
#'   that layer.
#' @param bw_mult,bw_floor,grid_n kernel density settings.
#' @param sigma_quantile quantile of the within-layer distances used as \eqn{\sigma_h}.
#' @param gamma exponent of the occupancy overlap (0 switches it off).
#' @return a numeric matrix of similarities (higher = more similar), `NULL` when no
#'   population passes the cell filter.
#' @examples
#' f <- system.file("extdata", "colon_fibroblast_two_sections.csv", package = "TICR")
#' d <- read.csv(f)
#' S <- population_similarity(d, groups = c("Day0", "Day9"))
#' round(S["Fibro 2", ], 2)
#' @export
population_similarity <- function(data, groups = NULL, layer = "layer", population = "population",
                                  r = "r", group = "condition", min_cells = 30,
                                  bw_mult = 2, bw_floor = 0.05, grid_n = 512,
                                  sigma_quantile = 0.5, gamma = 1) {
  stopifnot(all(c(layer, population, r, group) %in% names(data)))
  if (is.null(groups)) groups <- unique(data[[group]])[1:2]
  A <- data[data[[group]] == groups[1], , drop = FALSE]
  B <- data[data[[group]] == groups[2], , drop = FALSE]
  layers <- sort(unique(c(A[[layer]], B[[layer]])))

  per_layer <- function(ln) {
    dA <- A[A[[layer]] == ln, , drop = FALSE]; dB <- B[B[[layer]] == ln, , drop = FALSE]
    tA <- table(dA[[population]]); tB <- table(dB[[population]])
    rs <- sort(names(tA)[tA >= min_cells]); cs <- sort(names(tB)[tB >= min_cells])
    if (!length(rs) || !length(cs)) return(NULL)
    dA <- dA[dA[[population]] %in% rs, ]; dB <- dB[dB[[population]] %in% cs, ]
    rA <- split(dA[[r]], dA[[population]]); rB <- split(dB[[r]], dB[[population]])
    bw <- stats::bw.nrd0(c(dA[[r]], dB[[r]])) * bw_mult
    if (!is.finite(bw) || bw <= 0) bw <- bw_floor
    D <- matrix(NA_real_, length(rs), length(cs), dimnames = list(rs, cs))
    for (i in rs) for (j in cs) D[i, j] <- js_distance(rA[[i]], rB[[j]], bw, grid_n)
    sig <- as.numeric(stats::quantile(D, sigma_quantile, na.rm = TRUE))
    if (!is.finite(sig) || sig <= 0) sig <- 1
    list(S = exp(-D / sig), nA = tA[rs], nB = tB[cs])
  }
  J <- lapply(layers, per_layer); names(J) <- layers
  J <- J[!vapply(J, is.null, logical(1))]
  if (!length(J)) return(NULL)
  ar <- sort(unique(unlist(lapply(J, function(x) rownames(x$S)))))
  ac <- sort(unique(unlist(lapply(J, function(x) colnames(x$S)))))
  Sn <- matrix(0, length(ar), length(ac), dimnames = list(ar, ac)); Sd <- Sn
  for (ly in names(J)) {
    S <- J[[ly]]$S
    for (i in rownames(S)) for (j in colnames(S)) {
      s <- S[i, j]; if (!is.finite(s)) next
      w <- sqrt(as.numeric(J[[ly]]$nA[i]) * as.numeric(J[[ly]]$nB[j]))
      if (!is.finite(w) || w <= 0) next
      Sn[i, j] <- Sn[i, j] + w * s; Sd[i, j] <- Sd[i, j] + w
    }
  }
  Ss <- Sn / Sd; Ss[Sd == 0] <- NA
  if (gamma != 0) {
    occ <- function(d) { t <- table(d[[population]], d[[layer]]); t / rowSums(t) }
    pA <- occ(A); pB <- occ(B)
    LO <- matrix(NA_real_, nrow(Ss), ncol(Ss), dimnames = dimnames(Ss))
    for (i in rownames(Ss)) for (j in colnames(Ss)) {
      lys <- union(colnames(pA), colnames(pB))
      pi <- vapply(lys, function(l) if (l %in% colnames(pA)) pA[i, l] else 0, numeric(1))
      pj <- vapply(lys, function(l) if (l %in% colnames(pB)) pB[j, l] else 0, numeric(1))
      LO[i, j] <- sum(sqrt(pi * pj))
    }
    Ss <- Ss * LO^gamma
  }
  Ss
}

#' Rank the candidate populations of one query population
#'
#' @param S similarity matrix from [population_similarity()].
#' @param query a row name of `S`.
#' @param target optional column name; when given, its rank (1 = best, ties take the
#'   best rank) and normalized rank (0 = best) are returned as attributes.
#' @return a data.frame of candidates sorted by decreasing similarity.
#' @examples
#' f <- system.file("extdata", "colon_fibroblast_two_sections.csv", package = "TICR")
#' S <- population_similarity(read.csv(f), groups = c("Day0", "Day9"))
#' rank_candidates(S, "Fibro 2", target = "IAF 4")
#' @export
rank_candidates <- function(S, query, target = NULL) {
  stopifnot(query %in% rownames(S))
  v <- S[query, ]; v <- v[is.finite(v)]
  out <- data.frame(candidate = names(v), similarity = as.numeric(v), stringsAsFactors = FALSE)
  out <- out[order(-out$similarity), ]; rownames(out) <- NULL
  if (!is.null(target) && target %in% names(v)) {
    rk <- rank(-v, ties.method = "min")[target]
    attr(out, "rank") <- unname(rk)
    attr(out, "norm_rank") <- unname((rk - 1) / max(1, length(v) - 1))
  }
  out
}

#' Average similarity matrices over section pairs
#'
#' Element-wise mean over a list of matrices of possibly different dimensions; entries
#' present in fewer than `min_pairs` matrices are set to `NA`.
#'
#' @param mats list of matrices from [population_similarity()].
#' @param min_pairs minimum number of contributing pairs per entry.
#' @return a matrix with the union of all row and column names.
#' @export
pool_similarity <- function(mats, min_pairs = 1) {
  mats <- mats[!vapply(mats, is.null, logical(1))]
  rows <- sort(unique(unlist(lapply(mats, rownames))))
  cols <- sort(unique(unlist(lapply(mats, colnames))))
  acc <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols)); cnt <- acc
  for (S in mats) for (i in rownames(S)) for (j in colnames(S)) {
    v <- S[i, j]; if (is.finite(v)) { acc[i, j] <- acc[i, j] + v; cnt[i, j] <- cnt[i, j] + 1 }
  }
  M <- acc / cnt; M[cnt < min_pairs] <- NA; M
}
