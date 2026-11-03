# =====================================================================
# Internal MorphoGAM curve/coordinate fitter (ported verbatim from the colon pipeline; see ?fit_layer_coordinates)
# =====================================================================
# Design (settled after extensive testing; superseded methods deleted):
#
#   1. CONNECTIVITY GUARD (auto-k): raise the kNN k until the graph is one
#      connected component (bridge any leftover). This is the single change that
#      lets MorphoGAM run automatically on every slice - the original pipeline's
#      default knn=5 disconnected on sparse SM/ME layers.
#
#   2. OPEN vs CLOSED = MorphoGAM's OWN detector. MorphoGAM::CurveFinder(..., loop="auto")
#      uses a Mapper-graph cycle test (PC1 lens); fed our connected graph it is
#      reliable. We do NOT reinvent this. (Geometric detectors we tried all
#      failed on the non-convex/ruffled shapes; a phi1-lens Mapper was no better.)
#
#   3. H1 PERSISTENCE (ripserr) as an independent CROSS-CHECK. When it disagrees
#      with loop="auto", we flag the layer for review - those are the genuinely
#      ambiguous (thick/fat/fragmented) shapes, and loop mostly affects t, not the
#      r_for_jsd the benchmark uses.
#
#   4. knot.fx REGULARIZATION + BINNED FALLBACK. Fit with default knot.fx (detail);
#      if the curve spikes in a sparse gap, refit with knot.fx=40; if it still
#      spikes or errors, use a robust binned-centerline fit.
#
# Output is drop-in for the three-step engine: r_for_jsd in [0,1] + metadata.
# =====================================================================


# ---- REML-smoothed curve f(t): mgcv auto-selects wiggliness (REML+select), more
#      stable than MorphoGAM's default GCV (no pathological spikes). t is taken from
#      MorphoGAM's geodesic embedding; knn (not knot.fx) is the real lever. ----
.reml_curve <- function(t, xy, closed, kmax = 80) {
  k  <- min(kmax, max(10, round(0.3 * length(t))))
  bs <- if (closed) "cc" else "cr"
  fx <- mgcv::gam(xy[,1] ~ s(t, bs = bs, k = k), method = "REML", select = TRUE)
  fy <- mgcv::gam(xy[,2] ~ s(t, bs = bs, k = k), method = "REML", select = TRUE)
  cx <- as.numeric(predict(fx)); cy <- as.numeric(predict(fy))
  r  <- sqrt((xy[,1]-cx)^2 + (xy[,2]-cy)^2)                 # orthogonal dist to curve
  gc <- colMeans(xy); dcen <- sqrt((xy[,1]-gc[1])^2 + (xy[,2]-gc[2])^2)
  dctr <- sqrt((cx-gc[1])^2 + (cy-gc[2])^2)
  rs <- sign(dcen - dctr) * r                              # outer +, inner -
  if (suppressWarnings(cor(rs, dcen)) < 0) rs <- -rs
  list(cx = cx, cy = cy, r01 = (rs - min(rs))/(max(rs) - min(rs)),
       edf = round(sum(fx$edf) + sum(fy$edf), 1))
}

# ---- 1. connectivity guard ----
.knn_adj <- function(xy, k) {
  n <- nrow(xy); k <- min(k, n - 1)
  nn <- FNN::get.knn(xy, k = k)$nn.index
  i <- rep(seq_len(n), each = k); j <- as.vector(t(nn))
  A <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
  A <- A + Matrix::t(A); A@x <- rep(1, length(A@x)); A
}
.connect <- function(xy, k_init = 15, k_max = 100, k_step = 10) {
  k <- k_init
  repeat {
    A <- .knn_adj(xy, k)
    comp <- igraph::components(igraph::graph_from_adjacency_matrix(A, "undirected", diag = FALSE))
    if (comp$no == 1 || k >= k_max) break
    k <- k + k_step
  }
  if (comp$no > 1) {                                   # bridge leftover components
    memb <- comp$membership; ord <- order(comp$csize, decreasing = TRUE)
    base_idx <- which(memb == ord[1]); ei <- integer(0); ej <- integer(0)
    for (ci in ord[-1]) {
      pts <- which(memb == ci)
      nn  <- FNN::get.knnx(xy[base_idx,,drop=FALSE], xy[pts,,drop=FALSE], k = 1)
      w   <- which.min(nn$nn.dist[,1]); ei <- c(ei, pts[w]); ej <- c(ej, base_idx[nn$nn.index[w,1]])
      base_idx <- c(base_idx, pts)
    }
    B <- Matrix::sparseMatrix(i = c(ei,ej), j = c(ej,ei), x = 1, dims = dim(A))
    A <- A + B; A@x <- rep(1, length(A@x))
  }
  list(A = A, k = k)
}

# ---- 3. H1-persistence cross-check (closed if a loop persists across scales) ----
detect_closed_tda <- function(xy, subsample = 500, persist_thresh = 0.2) {
  if (!requireNamespace("ripserr", quietly = TRUE)) return(list(closed = NA, h1 = NA_real_))
  n <- nrow(xy); set.seed(1); idx <- if (n > subsample) sample(n, subsample) else seq_len(n)
  s <- xy[idx,,drop=FALSE]; diam <- sqrt(diff(range(s[,1]))^2 + diff(range(s[,2]))^2)
  ph <- tryCatch(as.data.frame(ripserr::vietoris_rips(s, max_dim = 1)), error = function(e) NULL)
  if (is.null(ph)) return(list(closed = FALSE, h1 = 0))
  h1 <- ph[ph$dimension == 1, ]; ps <- if (nrow(h1)) max(h1$death - h1$birth) / diam else 0
  list(closed = ps > persist_thresh, h1 = round(ps, 3))
}

# ---- spike = EVALUATION score (NOT a reject flag): how far the centerline strays into
#   EMPTY space. Considers ONLY the centerline points NOT on the cell cloud (nearest cell
#   > 1 cell-spacing away) - on-tissue points sit at ~0 and would just dilute the score.
#   Reports the mean of those off-cloud distances (in cell-spacing units) after dropping the
#   worst 5% (inevitable overshoots). 0 = whole curve sits on the cloud; higher = the curve
#   spends more of its length, deeper, in the void. Used to rank fits / slices, not to flag. ----
.SPIKE_THRESH <- 5     # retained only as a param default; spike no longer gates/flags
.spike_score <- function(xy, ctr) {
  cs <- median(FNN::get.knn(xy, k = 1)$nn.dist)
  d  <- FNN::get.knnx(xy, as.matrix(ctr), k = 1)$nn.dist[,1] / max(cs, 1e-9)   # spacing units
  off <- d[d > 1]                                  # points NOT on the cell cloud (>1 spacing)
  if (!length(off)) return(0)                      # curve sits entirely on the cloud
  mean(off[off <= stats::quantile(off, 0.95)])     # typical off-cloud excursion depth, trimmed
}

# ---- shape quality: is this a thin 1-D curve, or a filled/diffuse 2-D region
#      where the 1-D coordinate (r_for_jsd) is inherently unreliable? ----
#   inner_occ = fraction of cells in the inner 40% radius (filled center -> high)
#   thinness  = median per-angle radial IQR / median radius (thick/diffuse -> high)
.shape_quality <- function(xy) {
  gc <- colMeans(xy); rad <- sqrt((xy[,1]-gc[1])^2 + (xy[,2]-gc[2])^2)
  rmax <- as.numeric(quantile(rad, 0.98)); inner_occ <- mean(rad < 0.4 * rmax)
  ang <- atan2(xy[,2]-gc[2], xy[,1]-gc[1]); bin <- cut(ang, 36, labels = FALSE)
  spread <- tapply(rad, bin, function(r) diff(quantile(r, c(.25,.75))))
  thinness <- as.numeric(median(spread, na.rm = TRUE) / median(rad))
  list(inner_occ = round(inner_occ,3), thinness = round(thinness,2),
       diffuse = (inner_occ > 0.10) || (thinness > 0.25))
}

# ---- package a MorphoGAM fit into the three-step schema (r_for_jsd, lumen=0) ----
.package_mg <- function(fit, xy, method, loop) {
  r <- as.numeric(fit$xyt$r); gc <- colMeans(xy); dcen <- sqrt((xy[,1]-gc[1])^2 + (xy[,2]-gc[2])^2)
  if (suppressWarnings(cor(r, dcen)) < 0) r <- -r
  r01 <- (r - min(r)) / (max(r) - min(r))
  out <- data.frame(x = xy[,1], y = xy[,2], t = as.numeric(fit$xyt$t),
                    r = r01, r_original = r01, r_for_jsd = r01,
                    ctr_x = as.numeric(fit$xyt$f1), ctr_y = as.numeric(fit$xyt$f2))
  attr(out, "mode") <- if (isTRUE(loop)) "closed" else "open"; attr(out, "method") <- method; out
}

# ---- 4. robust binned-centerline fallback (given closed/open) ----
.binned_fit <- function(xy, closed, n_bins = 200, smooth_win = 3) {
  cg <- .connect(xy); d <- pmax(Matrix::rowSums(cg$A), 1e-9)
  Di <- Matrix::Diagonal(x = 1/sqrt(d)); S <- Di %*% cg$A %*% Di
  e <- RSpectra::eigs_sym(S, k = 3, which = "LA")
  phi1 <- e$vectors[,2]; phi2 <- e$vectors[,3]
  if (closed) { s <- atan2(phi2, phi1); t01 <- (s + pi)/(2*pi); breaks <- seq(-pi, pi, length.out = n_bins+1); periodic <- TRUE }
  else        { s <- phi1; t01 <- (rank(s)-1)/(length(s)-1); breaks <- seq(min(s), max(s), length.out = n_bins+1); periodic <- FALSE }
  bin <- cut(s, breaks = breaks, include.lowest = TRUE, labels = FALSE); B <- n_bins
  cx <- tapply(xy[,1], bin, mean); cy <- tapply(xy[,2], bin, mean)
  fx <- rep(NA_real_, B); fy <- rep(NA_real_, B); fx[as.integer(names(cx))] <- cx; fy[as.integer(names(cy))] <- cy
  fill <- function(v){ idx <- which(!is.na(v)); if(!length(idx)) return(v)
    for (b in which(is.na(v))){ dd <- if(periodic) pmin(abs(idx-b), B-abs(idx-b)) else abs(idx-b); v[b] <- v[idx[which.min(dd)]] }; v }
  sm <- function(v){ n<-length(v); out<-numeric(n); h<-smooth_win
    for (b in seq_len(n)){ w <- if(periodic) ((b-h-1):(b+h-1))%%n+1 else max(1,b-h):min(n,b+h); out[b]<-mean(v[w]) }; out }
  fx <- sm(fill(fx)); fy <- sm(fill(fy)); px <- fx; py <- fy
  ix1 <- if(periodic) c(2:B,1) else c(2:B,B); ix0 <- if(periodic) c(B,1:(B-1)) else c(1,1:(B-1))
  tx <- px[ix1]-px[ix0]; ty <- py[ix1]-py[ix0]; L <- sqrt(tx^2+ty^2); L[L==0]<-1; tx<-tx/L; ty<-ty/L
  nx <- -ty; ny <- tx
  rs <- (xy[,1]-px[bin])*nx[bin] + (xy[,2]-py[bin])*ny[bin]
  gc <- colMeans(xy); dcen <- sqrt((xy[,1]-gc[1])^2+(xy[,2]-gc[2])^2)
  if (suppressWarnings(cor(rs, dcen)) < 0) rs <- -rs
  r01 <- (rs - min(rs))/(max(rs)-min(rs))
  data.frame(x = xy[,1], y = xy[,2], t = t01, r = r01, r_original = r01, r_for_jsd = r01,
             ctr_x = px[bin], ctr_y = py[bin])
}

# ---- project ALL cells onto a fitted centerline polyline (cx,cy ordered by t).
#      Returns per-cell t, nearest centerline point, and SIGNED perpendicular
#      distance (local normal). Used so we can FIT the curve on a subsample (which
#      separates dense touching folds -> no geodesic short-circuit) yet still assign
#      a coordinate to every cell. ----
.project_all <- function(xy_all, cx, cy, tvec, closed) {
  ord <- order(tvec); px <- cx[ord]; py <- cy[ord]; pt <- tvec[ord]
  keep <- c(TRUE, (diff(px) != 0) | (diff(py) != 0))   # drop repeated vertices (stable tangents)
  px <- px[keep]; py <- py[keep]; pt <- pt[keep]; nP <- length(px)
  pt <- (pt - min(pt)) / (max(pt) - min(pt) + 1e-12)    # t -> [0,1]
  nn <- FNN::get.knnx(cbind(px, py), xy_all, k = 1)$nn.index[,1]
  if (isTRUE(closed)) { i1 <- nn %% nP + 1; i0 <- (nn - 2) %% nP + 1 }
  else                { i1 <- pmin(nn + 1, nP); i0 <- pmax(nn - 1, 1) }
  tx <- px[i1] - px[i0]; ty <- py[i1] - py[i0]; L <- sqrt(tx^2 + ty^2); L[L == 0] <- 1
  nx <- -ty / L; ny <- tx / L
  rs <- (xy_all[,1] - px[nn]) * nx + (xy_all[,2] - py[nn]) * ny
  list(t = pt[nn], ctr_x = px[nn], ctr_y = py[nn], rs = rs)
}

# ---- orient signed distance so the LUMINAL side is small: r=0 nearest the lumen
#      (whole-tissue `center`, passed by the caller), r=1 at the outer/serosal edge.
#      Using the global tissue center - not the per-layer centroid - makes the
#      inner/outer sign CONSISTENT across ring, arc and folded layers. ----
.orient01 <- function(rs, xy, center) {
  d <- sqrt((xy[,1] - center[1])^2 + (xy[,2] - center[2])^2)
  if (isTRUE(suppressWarnings(cor(rs, d)) < 0)) rs <- -rs   # far from lumen -> larger r
  (rs - min(rs)) / (max(rs) - min(rs) + 1e-12)              # inner (near lumen) -> 0
}

# ---- centerline polyline (x,y) ordered along t; close the ring if closed ----
.line_from <- function(cx, cy, t, closed) {
  o <- order(t); px <- cx[o]; py <- cy[o]
  keep <- c(TRUE, (diff(px) != 0) | (diff(py) != 0)); px <- px[keep]; py <- py[keep]
  if (isTRUE(closed) && length(px) > 2) { px <- c(px, px[1]); py <- c(py, py[1]) }
  data.frame(x = px, y = py)
}
# ---- line-shape detectors: CROSSING (self-intersection) & ZIGZAG (sharp reversal).
#      These catch bad centerlines that the spike score MISSES - a chord through
#      populated tissue keeps every curve point near a cell (low spike) yet crosses. ----
.resample_line <- function(line, m = 160) {
  n <- nrow(line); if (n <= m) return(line); line[round(seq(1, n, length.out = m)), , drop = FALSE] }
.seg_int <- function(p1,p2,p3,p4) {
  d <- function(a,b,c) (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1])
  d1<-d(p3,p4,p1); d2<-d(p3,p4,p2); d3<-d(p1,p2,p3); d4<-d(p1,p2,p4)
  ((d1>0&d2<0)|(d1<0&d2>0)) & ((d3>0&d4<0)|(d3<0&d4>0)) }
.line_selfcross <- function(line, m = 160) {
  L <- .resample_line(line, m); P <- as.matrix(L[,c("x","y")]); n <- nrow(P)-1
  if (n < 3) return(0L); cnt <- 0L
  for (i in 1:(n-2)) for (j in (i+2):n) { if (i==1 && j==n) next
    if (isTRUE(.seg_int(P[i,],P[i+1,],P[j,],P[j+1,]))) cnt <- cnt + 1L }; cnt }
.line_maxturn <- function(line, m = 160) {
  L <- .resample_line(line, m); P <- as.matrix(L[,c("x","y")]); n <- nrow(P); if (n < 3) return(0)
  v1 <- P[2:(n-1),,drop=FALSE]-P[1:(n-2),,drop=FALSE]; v2 <- P[3:n,,drop=FALSE]-P[2:(n-1),,drop=FALSE]
  cang <- pmin(1, pmax(-1, rowSums(v1*v2)/(sqrt(rowSums(v1^2))*sqrt(rowSums(v2^2))+1e-12)))
  round(max(acos(cang)*180/pi, na.rm=TRUE), 0) }
# verdict for a centerline: clean / CROSSING / ZIGZAG. Spike is NOT a flag (a spike is
# often inevitable) - it is reported as an evaluation score instead. `spike` arg unused.
.line_flag <- function(line, spike = NA) {
  if (.line_selfcross(line) > 0) return("CROSSING")
  if (.line_maxturn(line) > 150) return("ZIGZAG")
  "clean" }

# ---- build one fit-output df (project all cells onto a centerline, orient r, tag attrs) ----
.mk_out <- function(xy, cx, cy, tvec, closed, center, method, mg_plots = NULL, spike_val = NULL) {
  pr <- .project_all(xy, cx, cy, tvec, closed); r01 <- .orient01(pr$rs, xy, center)
  out <- data.frame(x = xy[,1], y = xy[,2], t = pr$t, r = r01,
                    r_original = r01, r_for_jsd = r01, ctr_x = pr$ctr_x, ctr_y = pr$ctr_y)
  ss <- if (is.null(spike_val)) .spike_score(xy, cbind(pr$ctr_x, pr$ctr_y)) else spike_val
  attr(out,"mode") <- if (closed) "closed" else "open"; attr(out,"method") <- method
  attr(out,"spike") <- round(ss,1)
  attr(out,"flag") <- .line_flag(.line_from(pr$ctr_x, pr$ctr_y, pr$t, closed), spike = round(ss,1))
  attr(out,"mg_plots") <- mg_plots; out
}

# fit_candidates - the SIMPLE-FIRST escalation, returning EVERY candidate it computed.
#   Step 1 plain MorphoGAM (default knn; crash -> auto-k). Step 2 keep plain if clean.
#   Step 3 (only if flagged) adaptive subsample+REML; then binned if it still spikes.
#   Returns list(chosen, chosen_key, plain, reml, binned) - each candidate a fit df (or
#   NULL if not computed). `chosen` is whichever the ladder selected. Crossing/zigzag are
#   flags (kept, not rejected past REML); only spike>thresh forces binned.
fit_candidates <- function(xy, spike_thresh = .SPIKE_THRESH, force_loop = NULL, center = NULL, fit_cap = 6000) {
  xy <- as.matrix(xy); if (is.null(center)) center <- colMeans(xy)
  q <- .shape_quality(xy); quality <- if (q$diffuse) "diffuse" else "clean"
  loop_arg <- if (is.null(force_loop)) "auto" else isTRUE(force_loop)
  loop_source <- if (is.null(force_loop)) "morphogam-auto" else "user-force"
  cand <- list(plain = NULL, reml = NULL, binned = NULL)
  finish <- function(o, key, kused) { attr(o,"loop_source")<-loop_source; attr(o,"quality")<-quality
    attr(o,"k_used")<-kused; attr(o,"cand_key")<-key; o }
  ret <- function(key) c(list(chosen = cand[[key]], chosen_key = key), cand)

  # ---- Step 1: plain MorphoGAM, default knn; crash -> auto-k ----
  fit <- tryCatch(MorphoGAM::CurveFinder(xy, loop = loop_arg), error = function(e) NULL); knn_lbl <- "default"
  if (is.null(fit)) { kk <- .connect(xy)$k
    fit <- tryCatch(MorphoGAM::CurveFinder(xy, knn = kk, loop = loop_arg), error = function(e) NULL)
    knn_lbl <- sprintf("auto-k=%d", kk) }
  if (!is.null(fit) && !is.null(fit$xyt)) {
    lp <- isTRUE(fit$loop); cx <- as.numeric(fit$xyt$f1); cy <- as.numeric(fit$xyt$f2); tt <- as.numeric(fit$xyt$t)
    p <- .mk_out(xy, cx, cy, tt, lp, center, sprintf("plain-mg(loop=%s,knn=%s)", ifelse(lp,"T","F"), knn_lbl),
                 mg_plots = list(curve = fit$curve.plot, coordinate = fit$coordinate.plot, residuals = fit$residuals.plot))
    cand$plain <- finish(p, "plain", knn_lbl)
    if (identical(attr(cand$plain,"flag"), "clean")) return(ret("plain"))   # Step 2: plain wins
  }

  # ---- Step 3: adaptive subsample + REML ----
  n <- nrow(xy); caps <- unique(pmin(n, c(fit_cap, round(fit_cap*0.66), round(fit_cap*0.42))))
  best <- NULL; best_key <- c(Inf, Inf)
  for (cap in caps) {
    set.seed(1); fi <- if (n > cap) sort(sample(n, cap)) else seq_len(n)
    xyf <- xy[fi, , drop = FALSE]; k <- .connect(xyf)$k
    fa <- tryCatch(MorphoGAM::CurveFinder(xyf, knn = k, loop = loop_arg), error = function(e) NULL)
    lpc <- if (is.null(fa)) (if (is.null(force_loop)) NA else isTRUE(force_loop)) else isTRUE(fa$loop)
    if (is.null(fa) || is.na(lpc)) next
    tf <- as.numeric(fa$xyt$t); rc <- .reml_curve(tf, xyf, lpc)
    pr <- .project_all(xy, rc$cx, rc$cy, tf, lpc); ss <- .spike_score(xy, cbind(pr$ctr_x, pr$ctr_y))
    cross <- .line_selfcross(.line_from(pr$ctr_x, pr$ctr_y, pr$t, lpc)); key <- c(cross, ss)
    if (is.null(best) || key[1] < best_key[1] || (key[1]==best_key[1] && key[2] < best_key[2])) {
      best <- list(cx=rc$cx, cy=rc$cy, tf=tf, lp=lpc, edf=rc$edf, cap=nrow(xyf), ss=ss, k=k); best_key <- key }
    if (cross == 0) break                                    # stop thinning once crossings clear
  }
  if (!is.null(best)) {
    r <- .mk_out(xy, best$cx, best$cy, best$tf, best$lp, center,
                 sprintf("resort:mg+REML(fitN=%d,edf=%.0f)", best$cap, best$edf), spike_val = best$ss)
    cand$reml <- finish(r, "reml", best$k)
    if (identical(attr(cand$reml, "flag"), "clean")) return(ret("reml"))   # REML clean -> keep it
    # else the REML fit STILL crosses/zigzags -> fall to binned (last-last resort)
  }
  # ---- binned = LAST-RESORT: reached only when REML still crosses/zigzags (or MorphoGAM+REML
  #      couldn't fit). Spike is NOT a trigger - it is an evaluation metric only. ----
  lp <- if (!is.null(best)) best$lp else if (is.null(force_loop)) FALSE else isTRUE(force_loop)
  xyf <- xy[if (n > fit_cap) { set.seed(1); sort(sample(n, fit_cap)) } else seq_len(n), , drop = FALSE]
  k <- .connect(xyf)$k; nb <- max(120, min(350, round(nrow(xyf)/12)))
  bf <- .binned_fit(xyf, lp, n_bins = nb)
  b <- .mk_out(xy, bf$ctr_x, bf$ctr_y, bf$t, lp, center, sprintf("resort:binned(%s,nb=%d)", ifelse(lp,"closed","open"), nb))
  cand$binned <- finish(b, "binned", k); return(ret("binned"))
}

# FitCurve - returns just the chosen fit (backward-compatible). See fit_candidates.
FitCurve <- function(xy, spike_thresh = .SPIKE_THRESH, force_loop = NULL, center = NULL, fit_cap = 6000)
  fit_candidates(xy, spike_thresh, force_loop, center, fit_cap)$chosen

# ---- assemble a raw fit df (from .mk_out) + the raw slice rows into the 28-col
#      manual `*_xyt_mapped.csv` format EXACTLY (same names & order = drop-in). ----
.to_manual <- function(f, s, layer) {
  g   <- function(nm) if (nm %in% names(s)) s[[nm]] else NA
  sid <- as.character(g("Slice_ID")[1]); if (is.na(sid) || sid == "") sid <- NA_character_
  pref <- if (!is.na(sid) && grepl("_D0_", sid)) "Healthy_all" else if (!is.na(sid) && grepl("_D9_", sid)) "DSS9_all" else NA
  kraw <- attr(f, "k_used")
  knn_num <- if (is.numeric(kraw)) as.integer(kraw)
             else if (grepl("auto-k=", as.character(kraw))) as.integer(sub(".*auto-k=", "", as.character(kraw)))
             else if (identical(as.character(kraw), "default")) 5L else NA_integer_
  out <- data.frame(
    x = f$x, y = f$y, t = f$t, r = f$r, f1 = f$ctr_x, f2 = f$ctr_y,
    cell_barcode = paste0(g("Slice_ID"), "_", g("Cell_ID")), layer = layer,
    Tier1 = g("Tier1"), Tier2 = g("Tier2"), Tier3 = g("Tier3"), Leiden_neigh = g("Leiden_neigh"),
    Sample_type = g("Sample_type"), Slice_ID = g("Slice_ID"), Cell_ID = g("Cell_ID"), Mouse_ID = g("Mouse_ID"),
    Technical_repeat_number = g("Technical_repeat_number"),
    source_file = if (!is.na(sid)) paste0(sid, ".csv") else NA,
    r_original = f$r_original, r_for_jsd = f$r_for_jsd,
    morphogam_sample_label = pref, morphogam_slice_id = sid, morphogam_layer = layer,
    morphogam_knn = knn_num, morphogam_manual = FALSE,
    morphogam_mode = attr(f, "mode"), morphogam_r_flip = FALSE,
    flip_updated_from_existing_csv = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
  for (a in c("mode","method","loop_source","quality","spike","flag","k_used","mg_plots","cand_key")) attr(out, a) <- attr(f, a)
  out
}

# fit_layer_xyt - fit one slice x layer, return the CHOSEN fit in 28-col manual format.
#   @param return_all if TRUE, also attaches attr "candidates" (named list of EVERY
#     candidate the ladder computed - plain/reml/binned - each in 28-col format with its
#     own attrs incl mg_plots) and attr "chosen_key". Lets the driver keep all versions
#     of a layer that climbed the ladder.
fit_layer_xyt <- function(slice_df, layer, min_cells = 30, force_loop = NULL, return_all = FALSE) {
  lay <- ifelse(grepl("^MU", slice_df$Leiden_neigh), "MU",
         ifelse(grepl("^SM", slice_df$Leiden_neigh), "SM",
         ifelse(grepl("^ME", slice_df$Leiden_neigh), "ME", NA_character_)))
  s <- slice_df[!is.na(lay) & lay == layer & is.finite(slice_df$x) & is.finite(slice_df$y), , drop = FALSE]
  keep <- names(which(table(s$Tier2) > min_cells)); s <- s[s$Tier2 %in% keep, , drop = FALSE]
  if (nrow(s) == 0) return(NULL)
  gx <- slice_df$x[is.finite(slice_df$x)]; gy <- slice_df$y[is.finite(slice_df$y)]
  gcen <- c(stats::median(gx), stats::median(gy))
  cc  <- fit_candidates(as.matrix(s[, c("x","y")]), force_loop = force_loop, center = gcen)
  out <- .to_manual(cc$chosen, s, layer)
  if (return_all) {
    alts <- list()
    for (nm in c("plain","reml","binned")) if (!is.null(cc[[nm]])) alts[[nm]] <- .to_manual(cc[[nm]], s, layer)
    attr(out, "candidates") <- alts; attr(out, "chosen_key") <- cc$chosen_key
  }
  out
}
