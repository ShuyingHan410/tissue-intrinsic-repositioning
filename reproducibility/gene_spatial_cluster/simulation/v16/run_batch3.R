## ============================================================
## Batch 3 — robustness: Experiments 5, 7, 7b, 8, 9, 10, 11, 6b.  Plan v16 §0.4.
## Reference design unless a cell says otherwise: k = 10, M = 64, NB size 1 (Fibro 6: size 2),
## matched pseudo-arms redrawn per replicate, gate SE0 < 0.05 & cells >= 30, R = 200.
## One target gene + 49 null backdrop genes per replicate (DESpace2 panel), all 5 methods,
## except Exp 9 (proposed with lambda, GAM-fs, DESpace2 + 300 blanks) and Exp 10 (proposed, GAM-fs).
##
## Usage: Rscript run_batch3.R <cell_id>   (1..N, see CELLS); no arg = all sequentially
## Env: R_N, R_FDR, CORES, SUBSTRATE_ROOT, BATCH3_OUT
## ============================================================
suppressMessages({source("engine_v16.R"); source("nd_generator.R"); source(file.path(PILOT, "null3_resampler.R"))
  library(DESpace); library(SingleCellExperiment); library(parallel)})
R_N <- as.integer(Sys.getenv("R_N", "200")); R_FDR <- as.integer(Sys.getenv("R_FDR", "50"))
CORES <- as.integer(Sys.getenv("CORES", Sys.getenv("SLURM_CPUS_PER_TASK", max(1L, detectCores() - 1L))))
ROOT <- Sys.getenv("SUBSTRATE_ROOT", "../../.."); OUTDIR <- Sys.getenv("BATCH3_OUT", "batch3_out"); dir.create(OUTDIR, showWarnings = FALSE)
M_MED <- 64; M_BLK <- 16; G_BG <- 50; B_RAW <- 300

x <- load_substrate(ROOT)
pool2 <- substrate_pool(x, "Fibro 2"); pool6 <- substrate_pool(x, "Fibro 6")
M6 <- 43   # Fibro 6 reference abundance: same all-panel measurement gives Fibro 2 = 42, Fibro 6 = 28 -> 64 * 28/42 (see BATCH3 notes)
d_real <- make_real_arms(x); r0r <- d_real$r[d_real$C == "0"]; r9r <- d_real$r[d_real$C == "1"]
sm <- readRDS(file.path(PILOT, "null3A_smoke.rds")); TL <- sm$TL     # frozen moment-balanced tilt (pop +0.20)
real_boot <- function() {                 # within-slice bootstrap of the 10 real Day-0 slices, arms by mouse (2 v 2)
  real <- x %>% filter(Tier2 == "Fibro 2", day == "Day 0") %>% group_by(Slice_ID) %>% filter(n() >= MIN_CELLS) %>% ungroup()
  mice <- unique(real$mouse); arm_of <- setNames(c(0, 0, 1, 1)[sample(4)], mice)
  d <- real %>% group_by(Slice_ID, mouse) %>% slice_sample(prop = 1, replace = TRUE) %>% ungroup() %>%
    transmute(slice = Slice_ID, mouse, C = arm_of[mouse], r)
  d$slice <- factor(d$slice); d$mouse <- factor(d$mouse); d$C <- factor(d$C); d }

## ---------------------------------------------------------------- cells
TAUS4 <- c(0, .05, .15, .30)
CELLS <- bind_rows(
  ## Exp 5: real geometry (D0 Fibro 2 vs D9 IAF 4) + bootstrap
  data.frame(exp = "5", sub = "tilt", scen = c("null3A", "null3B"), D = 0, tau = .15, M_B = M_MED, size = 1, k = 10, shape = "exp", sigma = 0, mkind = "normal"),
  data.frame(exp = "5", sub = "real", scen = c("gained", "reversed"), D = .05, tau = .15, M_B = M_MED, size = 1, k = 10, shape = "exp", sigma = 0, mkind = "normal"),
  expand.grid(exp = "5", sub = "boot", scen = c("null1", "gained"), D = NA, tau = c(0, .15), M_B = M_MED, size = 1, k = 5, shape = "exp", sigma = 0, mkind = "normal", stringsAsFactors = FALSE) %>% mutate(D = ifelse(scen == "gained", .05, 0)),
  ## Exp 7: abundance strata (M = 27, 331); primary tau .05, secondary 0/.15/.30
  expand.grid(exp = "7", sub = "pseudo", scen = c("null1", "gained"), D = NA, tau = TAUS4, M_B = c(27, 331), size = 1, k = 10, shape = "exp", sigma = 0, mkind = "normal", stringsAsFactors = FALSE) %>% mutate(D = ifelse(scen == "gained", .05, 0)),
  data.frame(exp = "7", sub = "pseudo", scen = "gained", D = .025, tau = .05, M_B = c(27, 331), size = 1, k = 10, shape = "exp", sigma = 0, mkind = "normal"),
  ## Exp 7b: per-arm abundance M_B/M_A in {1/3, 1, 3}, gained .05 + matching null, pseudo + real
  expand.grid(exp = "7b", sub = c("pseudo", "real"), scen = c("null1", "gained"), D = NA, tau = c(.05, .15), M_B = c(21, 64, 192), size = 1, k = 10, shape = "exp", sigma = 0, mkind = "normal", stringsAsFactors = FALSE) %>% mutate(D = ifelse(scen == "gained", .05, 0)),
  ## Exp 8: k = 20 (Null 1 exists from 1b; gained .05 across tau)
  data.frame(exp = "8", sub = "pseudo", scen = "gained", D = .05, tau = TAUS4, M_B = M_MED, size = 1, k = 20, shape = "exp", sigma = 0, mkind = "normal"),
  ## Exp 9: deployed pipeline, lambda on, sigma_mouse x tau x mouse-effect kind
  expand.grid(exp = "9", sub = "pseudo", scen = c("null1", "gained"), D = NA, tau = c(.05, .15), M_B = M_MED, size = 1, k = 10, shape = "exp", sigma = c(.07, .20), mkind = c("normal", "t3", "outlier"), stringsAsFactors = FALSE) %>% mutate(D = ifelse(scen == "gained", .05, 0)),
  ## Exp 10: Fibro 6 substrate, size 2, M = Fibro 6 median
  expand.grid(exp = "10", sub = "fibro6", scen = c("null1", "gained"), D = NA, tau = c(.05, .15), M_B = M6, size = 2, k = 10, shape = "exp", sigma = 0, mkind = "normal", stringsAsFactors = FALSE) %>% mutate(D = ifelse(scen == "gained", .05, 0)),
  ## Exp 11: awkward shapes at delta = .05
  data.frame(exp = "11", sub = "pseudo", scen = "gained", D = .05, tau = .15, M_B = M_MED, size = 1, k = 10, shape = c("local", "bimodal"), sigma = 0, mkind = "normal"),
  ## Exp 6b: FDR with correlated gene blocks
  data.frame(exp = "6b", sub = "pseudo", scen = "fdr", D = .05, tau = .15, M_B = M_MED, size = 1, k = 10, shape = "exp", sigma = 0, mkind = "normal"))
CELLS$id <- seq_len(nrow(CELLS))
cat("batch 3 cells:", nrow(CELLS), "\n")
if (nzchar(Sys.getenv("LIST"))) { options(width = 200); print(as.data.frame(CELLS[, c("id","exp","sub","scen","D","tau","M_B","size","k","shape","sigma","mkind")])); quit(save = "no") }

## ---------------------------------------------------------------- pieces
substrate_for <- function(sub, k) switch(sub,
  pseudo = design_k(pool2, k), fibro6 = design_k(pool6, k), real = d_real,
  tilt = make_pseudo_arms_null3(pool2, TL$w, mice_per_arm = 4, slices_per_mouse = c(2, 3, 3, 2)), boot = real_boot())
theta_for <- function(d, scen, D, sub) {
  rA <- d$r[d$C == "0"]; rB <- d$r[d$C == "1"]
  switch(scen,
    null1 = c(0, 0), null3B = c(0, 0),
    null3A = c(sm$th0A, sm$th9A),                                    # both arms at delta .05 relative to their own cells
    gained = c(0, solve_theta(if (sub %in% c("real", "tilt")) rB else rA, D)),
    reversed = c(solve_theta(rA, -D / 2), solve_theta(rB, D / 2)))
}
methods_on <- function(dd, cnt, d, which = c("proposed", "gam_naive", "gam_fs", "despace2", "pbulk_de")) {
  mp <- m_proposed(dd); out <- c(proposed = mp$pval)
  if ("gam_naive" %in% which) out["gam_naive"] <- m_gam(dd, "naive")
  if ("gam_fs" %in% which) out["gam_fs"] <- m_gam(dd, "fs")
  if ("despace2" %in% which) { ds <- NA_real_; invisible(capture.output(ds <- unname(m_despace2_panel(cnt, d)[1]))); out["despace2"] <- ds }
  if ("pbulk_de" %in% which) out["pbulk_de"] <- m_pseudobulk_de(dd)
  c(out, I2 = mean(c(mp$I2_0, mp$I2_9), na.rm = TRUE), k0 = mp$k0, k9 = mp$k9)
}
## standard replicate (Exp 5, 7, 7b, 8, 10, 11)
one_rep_std <- function(i, cl) {
  d <- substrate_for(cl$sub, cl$k)
  MA <- if (cl$exp == "7b") M_MED else cl$M_B      # Exp 7 / 10 are abundance STRATA (both arms); only 7b differs per arm
  Mv <- c(MA, cl$M_B)
  if (cl$shape == "exp") { tc <- theta_for(d, cl$scen, cl$D, cl$sub); dd <- sim_design(d, tc, Mv, cl$size, tau_theta = cl$tau) }
  else { rA <- d$r[d$C == "0"]; par <- if (cl$shape == "local") solve_local(rA, cl$D) else solve_bimodal(rA, cl$D)
         dd <- sim_design_u(d, c(0, par), Mv, cl$size, u_s = setNames(rnorm(nlevels(d$slice), 0, cl$tau), levels(d$slice)), shape = cl$shape) }
  cnt <- sim_panel(d, c(0, 0), G_BG, Mv[1], cl$size, tau_theta = cl$tau); cnt[1, ] <- dd$Y
  methods_on(dd, cnt, d, if (cl$sub == "fibro6") c("proposed", "gam_fs") else c("proposed", "gam_naive", "gam_fs", "despace2", "pbulk_de"))
}
## Exp 9: gene + 300 blanks share the SAME mouse effects (plan §3 rule), lambda from the blanks, applied to the gene
one_rep_9 <- function(i, cl) {
  d <- substrate_for("pseudo", cl$k); b_m <- draw_mouse_effects(levels(d$mouse), cl$sigma, cl$mkind)
  us <- function() setNames(rnorm(nlevels(d$slice), 0, cl$tau), levels(d$slice))
  tc <- theta_for(d, cl$scen, cl$D, "pseudo")
  dd <- sim_design_u(d, tc, M_MED, cl$size, u_s = us(), b_m = b_m)
  bp <- vapply(seq_len(B_RAW), function(b) m_proposed(sim_design_u(d, c(0, 0), M_BLK, cl$size, u_s = us(), b_m = b_m))$pval, numeric(1))
  L <- lambda_of(bp); mp <- m_proposed(dd)
  cnt <- sim_panel(d, c(0, 0), G_BG, M_MED, cl$size, tau_theta = cl$tau); cnt[1, ] <- dd$Y
  ds <- NA_real_; invisible(capture.output(ds <- unname(m_despace2_panel(cnt, d)[1])))
  c(proposed = mp$pval, proposed_lambda = unname(apply_lambda(mp$pval, L[["applied"]])), gam_fs = m_gam(dd, "fs"), despace2 = ds,
    lam_raw = unname(L["raw"]), lam_applied = unname(L["applied"]), n_usable = unname(L["n"]), blank_fp = mean(bp < .05, na.rm = TRUE),
    blank_z_sd = sd(qnorm(bp[bp > 0 & bp < 1] / 2)), I2 = mean(c(mp$I2_0, mp$I2_9), na.rm = TRUE), k0 = mp$k0, k9 = mp$k9,
    blank_p_q05 = unname(quantile(bp, .05, na.rm = TRUE)), blank_p_q01 = unname(quantile(bp, .01, na.rm = TRUE)), blank_p_q001 = unname(quantile(bp, .001, na.rm = TRUE)))
}
## Exp 6b: FDR mixture with block-correlated slice nudges (rho = .7, blocks of 20)
one_rep_6b <- function(i, cl) {
  d <- substrate_for("pseudo", cl$k); G <- 500; G1 <- 100; truth <- c(rep(TRUE, G1), rep(FALSE, G - G1))
  TH1 <- c(0, solve_theta(pool2$r, .05)); U <- draw_block_nudges(G, levels(d$slice), cl$tau, rho = .7, block = 20)
  cnt <- matrix(0L, G, nrow(d), dimnames = list(paste0("g", 1:G), paste0("c", seq_len(nrow(d)))))
  for (g in 1:G) cnt[g, ] <- sim_design_u(d, if (truth[g]) TH1 else c(0, 0), M_MED, cl$size, u_s = U[g, ])$Y
  pv <- mclapply(1:G, function(g) { d$Y <- cnt[g, ]; c(proposed = m_proposed(d)$pval, gam_fs = m_gam(d, "fs"), pbulk_de = m_pseudobulk_de(d)) }, mc.cores = CORES, mc.set.seed = TRUE)
  P <- do.call(rbind, pv); ds <- NULL; invisible(capture.output(ds <- unname(m_despace2_panel(cnt, d)[rownames(cnt)]))); list(P = cbind(P, despace2 = ds), truth = truth)
}

run_cell <- function(id) {
  cl <- CELLS[id, ]; key <- sprintf("exp%s|%s|%s|D%.3f|tau%.2f|MB%d|k%d|%s|s%.2f|%s", cl$exp, cl$sub, cl$scen, cl$D, cl$tau, cl$M_B, cl$k, cl$shape, cl$sigma, cl$mkind)
  f <- file.path(OUTDIR, sprintf("cell_%03d.rds", id)); if (file.exists(f)) { cat("skip", key, "\n"); return(invisible(NULL)) }
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260920 + id); cat("START", key, " cores =", CORES, "\n"); flush.console(); t0 <- Sys.time()
  if (cl$exp == "6b") { reps <- lapply(seq_len(R_FDR), function(i) { r <- one_rep_6b(i, cl); cat(sprintf("  6b rep %d/%d\n", i, R_FDR)); flush.console(); r })
    out <- list(cell = cl, key = key, P = lapply(reps, `[[`, "P"), truth = reps[[1]]$truth, n_rep = R_FDR)
  } else {
    fn <- if (cl$exp == "9") one_rep_9 else one_rep_std
    res <- mclapply(seq_len(R_N), fn, cl = cl, mc.cores = CORES, mc.set.seed = TRUE)
    ok <- vapply(res, is.numeric, logical(1)); P <- do.call(rbind, res[ok]); out <- list(cell = cl, key = key, P = P, n_rep = sum(ok))
  }
  out$secs <- as.numeric(Sys.time() - t0, units = "secs"); saveRDS(out, f)
  if (cl$exp != "6b") { pc <- intersect(colnames(out$P), c("proposed", "proposed_lambda", "gam_naive", "gam_fs", "despace2", "pbulk_de")); rj <- apply(out$P[, pc, drop = FALSE] < .05, 2, mean, na.rm = TRUE)
    cat(sprintf("DONE  %s  %d reps  %.0f s  reject@0.05: %s\n", key, out$n_rep, out$secs, paste(sprintf("%s=%.3f", names(rj), rj), collapse = " "))) } else cat(sprintf("DONE  %s  %.0f s\n", key, out$secs))
  flush.console()
}
args <- commandArgs(trailingOnly = TRUE); for (id in if (length(args)) as.integer(args[1]) else CELLS$id) run_cell(id)
