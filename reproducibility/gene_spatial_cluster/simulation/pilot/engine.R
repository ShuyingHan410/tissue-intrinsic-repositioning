## =====================================================================
## Pilot engine for revised_gene_repositioning_simulation_plan_v10.md
## All functions; no side effects on load.
## =====================================================================
suppressMessages({
  library(readr); library(dplyr); library(purrr); library(tidyr)
  library(metafor); library(mgcv)
})

CHI1_MED <- qchisq(0.5, 1)          # 0.4549364
MIN_CELLS <- 30L
SE_MAX    <- 0.05

## ---------------------------------------------------------------- data
load_substrate <- function(root = "../../..") {
  sl <- suppressMessages(read_csv(file.path(root, "gene_spatial_cluster",
                                            "gene_spatial_lineage_slices.csv"),
                                  show_col_types = FALSE))
  map <- sl %>% dplyr::distinct(Slice_ID, mouse, day)
  fs <- list.files(file.path(root, "output_threelayer", "per_slice_morphogam"),
                   "_MU_xyt_mapped.csv$", full.names = TRUE)
  co <- map_dfr(fs, function(f)
    suppressMessages(read_csv(f, show_col_types = FALSE)) %>%
      dplyr::transmute(Slice_ID = as.character(Slice_ID), r = r_for_jsd, Tier2))
  co %>% dplyr::inner_join(map, by = "Slice_ID")
}

## Pooled Day-0 Fibro 2 coordinates + observed slice sizes
substrate_pool <- function(x, tier = "Fibro 2", dy = "Day 0") {
  z <- x %>% dplyr::filter(Tier2 == tier, day == dy)
  sizes <- z %>% dplyr::count(Slice_ID) %>% dplyr::filter(n >= MIN_CELLS) %>% dplyr::pull(n)
  list(r = z$r, sizes = sizes,
       n_slices = length(sizes),
       n_mice = n_distinct(z$mouse[z$Slice_ID %in%
                (z %>% dplyr::count(Slice_ID) %>% dplyr::filter(n >= MIN_CELLS) %>% dplyr::pull(Slice_ID))]))
}

## ------------------------------------------------- pseudo-arm substrate
## v10 §1: randomisation unit is the MOUSE. Synthetic pseudo-mice built by
## resampling the pooled Day-0 coordinate distribution.
make_pseudo_arms <- function(pool, mice_per_arm = 4,
                             slices_per_mouse = c(1, 1, 1, 2), cells = NULL) {
  stopifnot(length(slices_per_mouse) == mice_per_arm)
  out <- list(); k <- 0L
  for (arm in c(0, 1)) for (m in seq_len(mice_per_arm)) {
    for (j in seq_len(slices_per_mouse[m])) {
      k <- k + 1L
      n <- if (is.null(cells)) sample(pool$sizes, 1) else cells
      out[[k]] <- data.frame(
        slice = sprintf("A%d_M%d_S%d", arm, m, j),
        mouse = sprintf("A%d_M%d", arm, m),
        C     = arm,
        r     = sample(pool$r, n, replace = TRUE))
    }
  }
  d <- dplyr::bind_rows(out)
  d$slice <- factor(d$slice); d$mouse <- factor(d$mouse); d$C <- factor(d$C)
  d
}

## Real mismatched geometry (Day 0 Fibro 2 vs Day 9 IAF 4) -- for B2
make_real_arms <- function(x, cells_cap = NULL) {
  g <- function(tier, dy, arm) x %>% dplyr::filter(Tier2 == tier, day == dy) %>%
    dplyr::group_by(Slice_ID) %>% dplyr::filter(dplyr::n() >= MIN_CELLS) %>% dplyr::ungroup() %>%
    dplyr::transmute(slice = Slice_ID, mouse, C = arm, r)
  d <- dplyr::bind_rows(g("Fibro 2", "Day 0", 0), g("IAF 4", "Day 9", 1))
  if (!is.null(cells_cap))
    d <- d %>% dplyr::group_by(slice) %>% dplyr::slice_sample(n = min(cells_cap, n())) %>% dplyr::ungroup()
  d$slice <- factor(d$slice); d$mouse <- factor(d$mouse); d$C <- factor(d$C)
  d
}

## --------------------------------------------------- truth calibration
zscore <- function(r) (r - mean(r)) / sd(r)

## expected delta for weight h = exp(theta*z)
delta_exp <- function(r, theta) {
  z <- zscore(r); h <- exp(theta * z)
  sum(h * r) / sum(h) - mean(r)
}
solve_theta <- function(r, target, lo = -8, hi = 8) {
  if (abs(target) < 1e-12) return(0)
  f <- function(t) delta_exp(r, t) - target
  if (f(lo) * f(hi) > 0) return(NA_real_)
  uniroot(f, c(lo, hi), tol = 1e-10)$root
}

## ------------------------------------------------------------ simulate
## mu_i = M * h_i / sum(h_i)  -> total abundance is exactly M for any h
sim_slice <- function(r, theta, M, size = 2) {
  z <- zscore(r); h <- exp(theta * z)
  mu <- M * h / sum(h)
  rnbinom(length(r), mu = mu, size = size)
}

## Dirichlet-multinomial variant for Ablation A (fixes N_event, varies N_eff)
sim_slice_dm <- function(r, theta, N_event, alpha) {
  z <- zscore(r); h <- exp(theta * z); p <- h / sum(h)
  q <- rgamma(length(p), shape = alpha * p, rate = 1)
  if (sum(q) <= 0) q <- rep(1, length(p))
  q <- q / sum(q)
  as.vector(rmultinom(1, N_event, q))
}

## Simulate a whole design. theta_c is a length-2 vector (arm0, arm9).
## M may be a scalar (same abundance both arms) or length-2 c(M0, M9) --
## the per-arm form exists for the DE-only specificity scenario.
sim_design <- function(d, theta_c, M, size = 2,
                       sigma_mouse = 0, tau_theta = 0, dm = NULL) {
  Mv <- rep_len(M, 2)
  mice <- levels(d$mouse)
  bm <- setNames(rnorm(length(mice), 0, sigma_mouse), mice)
  sls <- levels(d$slice)
  us <- setNames(rnorm(length(sls), 0, tau_theta), sls)
  d$Y <- NA_integer_
  for (s in sls) {
    i <- which(d$slice == s)
    arm <- as.integer(as.character(d$C[i[1]])) + 1L
    th <- theta_c[arm] + bm[[as.character(d$mouse[i[1]])]] + us[[s]]
    d$Y[i] <- if (is.null(dm)) sim_slice(d$r[i], th, Mv[arm], size)
              else sim_slice_dm(d$r[i], th, dm$N_event, dm$alpha)
  }
  d
}

## ------------------------------------------------- the proposed method
## delta and its permutation SE0 per slice
slice_delta <- function(r, Y) {
  n <- length(r); tot <- sum(Y)
  if (tot <= 0) return(c(delta = NA, se = NA, Neff = NA, mol = 0))
  Neff <- tot^2 / sum(Y^2)
  d <- sum(Y * r) / tot - mean(r)
  v <- var(r) * (1 / Neff - 1 / n)
  c(delta = d, se = if (is.finite(v) && v > 0) sqrt(v) else NA, Neff = Neff, mol = tot)
}

per_slice_table <- function(d) {
  d %>% dplyr::group_by(slice, mouse, C) %>%
    dplyr::summarise(n = n(), as.data.frame(t(slice_delta(r, Y))), .groups = "drop") %>%
    dplyr::rename(delta = delta, se = se)
}

## REML + Knapp-Hartung pooling of one arm
pool_arm <- function(tb, kh = TRUE) {
  tb <- tb %>% dplyr::filter(is.finite(delta), is.finite(se), se > 0, se < SE_MAX, n >= MIN_CELLS)
  k <- nrow(tb)
  if (k < 2) return(list(mu = NA, se = NA, k = k, I2 = NA, p = NA))
  fit <- try(rma(yi = tb$delta, sei = tb$se, method = "REML",
                 test = if (kh) "knha" else "z"), silent = TRUE)
  if (inherits(fit, "try-error")) return(list(mu = NA, se = NA, k = k, I2 = NA, p = NA))
  list(mu = as.numeric(fit$beta), se = fit$se, k = k,
       I2 = fit$I2, p = fit$pval)
}

## fixed-effect inverse-variance pooling (Ablation B: tau^2 = 0)
pool_arm_fixed <- function(tb) {
  tb <- tb %>% dplyr::filter(is.finite(delta), is.finite(se), se > 0, se < SE_MAX, n >= MIN_CELLS)
  k <- nrow(tb)
  if (k < 2) return(list(mu = NA, se = NA, k = k, p = NA))
  w <- 1 / tb$se^2
  mu <- sum(w * tb$delta) / sum(w); se <- sqrt(1 / sum(w))
  list(mu = mu, se = se, k = k, p = 2 * pnorm(-abs(mu / se)))
}

## D1: contrast of the two arms, Welch-Satterthwaite df
d1_test <- function(a0, a9) {
  if (!is.finite(a0$se) || !is.finite(a9$se)) return(c(D = NA, se = NA, p = NA))
  D <- a9$mu - a0$mu; se <- sqrt(a0$se^2 + a9$se^2)
  nu <- se^4 / (a9$se^4 / max(a9$k - 1, 1) + a0$se^4 / max(a0$k - 1, 1))
  c(D = D, se = se, p = 2 * pt(-abs(D / se), df = max(nu, 1)))
}

## full proposed pipeline -> D1 p-value  (naive_se toggles Ablation A)
m_proposed <- function(d, kh = TRUE, fixed = FALSE, naive_se = FALSE) {
  tb <- per_slice_table(d)
  if (naive_se) {
    ns <- d %>% dplyr::group_by(slice) %>%
      dplyr::summarise(se_naive = sd(r) / sqrt(sum(Y)), .groups = "drop")
    tb <- tb %>% dplyr::left_join(ns, by = "slice") %>% dplyr::mutate(se = se_naive)
  }
  f <- if (fixed) pool_arm_fixed else function(z) pool_arm(z, kh)
  a0 <- f(tb %>% dplyr::filter(C == "0")); a9 <- f(tb %>% dplyr::filter(C == "1"))
  r <- d1_test(a0, a9)
  list(p = unname(r["D"] ) , pval = unname(r["p"]), D = unname(r["D"]),
       k0 = a0$k, k9 = a9$k, I2_0 = a0$I2 %||% NA, I2_9 = a9$I2 %||% NA, tb = tb)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## ------------------------------------------------------ GAM comparators
## NOTE: C must be an ORDERED factor here. With an unordered factor, s(r, by = C)
## fits one smooth PER LEVEL, so testing it asks "is there a gradient?" rather than
## "does the gradient differ between conditions?" -- which makes any stably
## positioned gene a guaranteed rejection. Null 2 caught this. With an ordered
## factor the s(r):Co term IS the difference smooth s_Delta(r).
m_gam <- function(d, form = c("naive", "fs")) {
  form <- match.arg(form)
  d$Co <- as.ordered(d$C)
  f <- if (form == "naive")
    Y ~ C + s(r, k = 8) + s(r, by = Co, k = 8)
  else
    Y ~ C + s(r, k = 8) + s(r, by = Co, k = 8) + s(r, slice, bs = "fs", k = 5) + s(mouse, bs = "re")
  m <- try(bam(f, family = nb(), data = d, method = "fREML", discrete = TRUE), silent = TRUE)
  if (inherits(m, "try-error")) return(NA_real_)
  st <- summary(m)$s.table
  i <- grep(":Co", rownames(st), fixed = TRUE)
  if (!length(i)) return(NA_real_)
  suppressWarnings(min(st[i, "p-value"], na.rm = TRUE))
}

## ---------------------------------------------- DESpace2 (REAL package)
## DESpace 2.4.0, dsp_test(). Operates on a GENE PANEL: edgeR needs cross-gene
## dispersion, so testing one gene at a time is not how the method works (and was
## the cause of the inflated Type I error in the first pilot). One call returns a
## p-value per gene. test="QLF" (quasi-likelihood F) is the package default.
m_despace2_panel <- function(cnt, d, K = 5, test = "QLF") {
  cd <- data.frame(
    slice = d$slice, C = factor(d$C),
    dom = factor(cut(d$r, seq(0, 1, length.out = K + 1),
                     include.lowest = TRUE, labels = FALSE)))
  rownames(cd) <- colnames(cnt)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = cnt), colData = cd)
  r <- try(suppressMessages(DESpace::dsp_test(
    sce, cluster_col = "dom", sample_col = "slice", condition_col = "C",
    test = test, verbose = FALSE)), silent = TRUE)
  if (inherits(r, "try-error")) return(setNames(rep(NA_real_, nrow(cnt)), rownames(cnt)))
  out <- setNames(rep(NA_real_, nrow(cnt)), rownames(cnt))
  out[as.character(r$gene_id)] <- r$PValue
  out
}

## Simulate a panel of G genes on one design, all with the same theta_c.
sim_panel <- function(d, theta_c, G, M, size = 2, sigma_mouse = 0, tau_theta = 0) {
  cnt <- matrix(0L, G, nrow(d),
                dimnames = list(paste0("g", seq_len(G)), paste0("c", seq_len(nrow(d)))))
  for (g in seq_len(G))
    cnt[g, ] <- sim_design(d, theta_c, M, size, sigma_mouse, tau_theta)$Y
  cnt
}

## ------------------------------------------- DESpace2-style (pseudobulk)
## X_gsd = sum of counts in slice s, domain d; NB/quasi-Poisson GLM
## with domain x condition interaction, slice as the replicate unit.
m_despace2 <- function(d, K = 5) {
  br <- seq(0, 1, length.out = K + 1)
  pb <- d %>% dplyr::mutate(dom = cut(r, br, include.lowest = TRUE, labels = FALSE)) %>%
    dplyr::group_by(slice, C, dom) %>% dplyr::summarise(X = sum(Y), ncell = n(), .groups = "drop") %>%
    dplyr::filter(ncell > 0)
  pb$dom <- factor(pb$dom)
  if (n_distinct(pb$dom) < 2) return(NA_real_)
  m1 <- try(glm(X ~ offset(log(ncell)) + dom * C, family = quasipoisson, data = pb), silent = TRUE)
  m0 <- try(glm(X ~ offset(log(ncell)) + dom + C, family = quasipoisson, data = pb), silent = TRUE)
  if (inherits(m1, "try-error") || inherits(m0, "try-error")) return(NA_real_)
  a <- try(anova(m0, m1, test = "F"), silent = TRUE)
  if (inherits(a, "try-error")) return(NA_real_)
  p <- a$`Pr(>F)`[2]; if (is.null(p) || !is.finite(p)) NA_real_ else p
}

## ------------------------------------------------ pseudobulk DE (abundance)
m_pseudobulk_de <- function(d) {
  pb <- d %>% dplyr::group_by(slice, C) %>% dplyr::summarise(X = sum(Y), ncell = n(), .groups = "drop")
  if (n_distinct(pb$C) < 2) return(NA_real_)
  m <- try(glm(X ~ offset(log(ncell)) + C, family = quasipoisson, data = pb), silent = TRUE)
  if (inherits(m, "try-error")) return(NA_real_)
  s <- summary(m)$coefficients
  i <- grep("^C", rownames(s)); if (!length(i)) return(NA_real_)
  s[i[1], 4]
}

## ------------------------------------------------------------- lambda
## v6 §3: floored, with the n<20 fallback
lambda_of <- function(p, min_n = 20) {
  p <- p[is.finite(p) & p > 0 & p < 1]
  if (length(p) < min_n) return(c(raw = NA, applied = 1, n = length(p)))
  raw <- median(qnorm(p / 2)^2) / CHI1_MED
  c(raw = raw, applied = max(1, raw), n = length(p))
}
apply_lambda <- function(p, lam) { z <- abs(qnorm(p / 2)); 2 * pnorm(-z / sqrt(lam)) }

## C2 tail calibration: order statistic k = ceil(0.95 (B+1))
tail_cut <- function(p_blank, q = 0.95) {
  z <- abs(qnorm(p_blank[is.finite(p_blank) & p_blank > 0 & p_blank < 1] / 2))
  if (!length(z)) return(NA_real_)
  B <- length(z); k <- min(B, ceiling(q * (B + 1)))
  sort(z)[k]
}

## ------------------------------------------------- comparator-truth audit
## Audit the GENERATIVE truth, not a realised sample: takes the two arms'
## coordinate DISTRIBUTIONS (r0, r9) and their thetas. On the matched
## substrate r0 and r9 are the same pool, so any difference comes from h alone.
audit_truth <- function(r0, r9, theta_c, K = 5, tol = 1e-6) {
  dt <- theta_c[2] - theta_c[1]
  zz <- range(c(zscore(r0), zscore(r9)))
  swing <- if (abs(dt) < tol) 0 else abs(exp(dt * zz[2]) / exp(dt * zz[1]) - 1)
  br <- seq(0, 1, length.out = K + 1)
  dom_p <- function(r, th) {
    w <- exp(th * zscore(r)); d <- cut(r, br, include.lowest = TRUE, labels = FALSE)
    s <- tapply(w, factor(d, levels = seq_len(K)), sum); s[is.na(s)] <- 0
    as.numeric(s / sum(s))
  }
  p0 <- dom_p(r0, theta_c[1]); p9 <- dom_p(r9, theta_c[2])
  maxdiff <- max(abs(p9 - p0))
  ## GAM null: h9/h0 constant in r. True iff theta equal AND same population.
  same_pop <- isTRUE(all.equal(sort(r0), sort(r9), tolerance = 1e-12))
  list(dtheta = dt, ratio_swing = swing,
       gam_null = (abs(dt) < tol) && same_pop,
       max_domain_diff = maxdiff,
       despace_null = maxdiff < 1e-4,
       domains = data.frame(dom = seq_len(K), arm0 = p0, arm1 = p9, diff = p9 - p0))
}

## Audit the REALISED substrate (the actual sampled coordinates), not the
## generative pool. Auditing r0 vs r9 when they are the same vector by
## construction returns 0 trivially and tests nothing -- the realised draw is
## what every method actually sees.
audit_realised <- function(d, K = 5) {
  br <- seq(0, 1, length.out = K + 1)
  sh <- function(a) {
    r <- d$r[d$C == a]
    as.numeric(table(cut(r, br, include.lowest = TRUE))) / length(r)
  }
  p0 <- sh("0"); p9 <- sh("1")
  list(max_domain_diff = max(abs(p9 - p0)),
       ks = suppressWarnings(ks.test(d$r[d$C == "0"], d$r[d$C == "1"])$statistic),
       domains = data.frame(dom = seq_len(K), arm0 = p0, arm1 = p9, diff = p9 - p0))
}

## ----------------------------------------------------------- MC helpers
## cluster-robust MC SE from per-replicate rejection proportions
mc_rate <- function(rej_by_rep) {
  pr <- vapply(rej_by_rep, function(v) mean(v, na.rm = TRUE), numeric(1))
  pr <- pr[is.finite(pr)]
  R <- length(pr); a <- mean(pr); se <- sd(pr) / sqrt(R)
  neff <- if (is.finite(se) && se > 0) a * (1 - a) / se^2 else NA
  c(alpha = a, se = se, lo = max(0, a - 1.96 * se), hi = min(1, a + 1.96 * se),
    R = R, n_eff = neff)
}
