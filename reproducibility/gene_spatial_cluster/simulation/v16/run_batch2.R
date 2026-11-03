## ============================================================
## Batch 2 — Experiments 2 (directional power), 3 (specificity), 6 (FDR).
## Plan v16 §0.4–0.5. Submit ONLY after summarise_exp1.R has been inspected.
## Reference design: k = 10, M = 64, NB size 1, matched pseudo-arms redrawn per
## replicate, gate SE0 < 0.05 & cells >= 30, R = 200 (FDR: R = 50).
## Constructions reproduce pilot/run_b1.R, run_spec.R, run_fdr.R exactly except
## size (1, not 2) and k (10; FDR also 20).
##
## Usage: Rscript run_batch2.R <cell_id>    cell_id in 1..30 (see CELLS)
##        Rscript run_batch2.R              all cells sequentially
## Env: R_N (200), R_FDR (50), G_FDR (500), G1_FDR (100), G_BG (50), CORES, SUBSTRATE_ROOT, BATCH2_OUT
## ============================================================
suppressMessages({source("engine_v16.R"); source("nd_generator.R"); library(DESpace); library(SingleCellExperiment); library(parallel)})
R_N   <- as.integer(Sys.getenv("R_N", "200"))
R_FDR <- as.integer(Sys.getenv("R_FDR", "50"))
G_FDR <- as.integer(Sys.getenv("G_FDR", "500")); G1_FDR <- as.integer(Sys.getenv("G1_FDR", "100"))
G_BG  <- as.integer(Sys.getenv("G_BG", "50"))
M_MED <- 64; SIZE <- 1
CORES <- as.integer(Sys.getenv("CORES", Sys.getenv("SLURM_CPUS_PER_TASK", max(1L, detectCores() - 1L))))
ROOT  <- Sys.getenv("SUBSTRATE_ROOT", "../../..")
OUTDIR <- Sys.getenv("BATCH2_OUT", "batch2_out"); dir.create(OUTDIR, showWarnings = FALSE)

TAUS <- c(0, 0.05, 0.15, 0.30); DS <- c(0.025, 0.05, 0.075)
CELLS <- rbind(
  data.frame(exp = "2", k = 10, scen = "gained",   D = rep(DS, each = 4), tau = rep(TAUS, 3), kappa = NA, Mratio = 1),
  data.frame(exp = "2", k = 10, scen = "reversed", D = rep(DS, each = 4), tau = rep(TAUS, 3), kappa = NA, Mratio = 1),
  data.frame(exp = "3", k = 10, scen = c("de_only", "repo_only", "nd_moderate", "nd_strong"),
             D = c(0, 0.05, 0, 0), tau = c(0.15, 0.15, 0, 0), kappa = c(NA, NA, 10, 50), Mratio = c(2, 1, 1, 1)),
  data.frame(exp = "6", k = c(10, 20), scen = "fdr", D = 0.05, tau = 0.15, kappa = NA, Mratio = 1))
CELLS$id <- seq_len(nrow(CELLS))

x <- load_substrate(ROOT); pool <- substrate_pool(x)
th <- function(target) solve_theta(pool$r, target)
theta_for <- function(scen, D) switch(scen,
  gained = , repo_only = c(0, th(D)),
  reversed = c(th(-D / 2), th(D / 2)),
  de_only = c(0, 0))

methods_on <- function(dd, cnt, d) {
  mp <- m_proposed(dd)
  ds <- NA_real_; invisible(capture.output(ds <- unname(m_despace2_panel(cnt, d)[1])))
  c(proposed = mp$pval, gam_naive = m_gam(dd, "naive"), gam_fs = m_gam(dd, "fs"), despace2 = ds,
    pbulk_de = m_pseudobulk_de(dd), I2 = mean(c(mp$I2_0, mp$I2_9), na.rm = TRUE),
    k0 = mp$k0, k9 = mp$k9)                                   # testability (§0b B1)
}
## Exp 2 / Exp 3 standard cells: one target gene + G_BG - 1 NULL backdrop genes
one_rep_std <- function(i, k, theta_c, Mv, tau) {
  d <- design_k(pool, k)
  dd <- sim_design(d, theta_c, Mv, SIZE, tau_theta = tau)
  cnt <- sim_panel(d, c(0, 0), G_BG, M_MED, SIZE, tau_theta = tau); cnt[1, ] <- dd$Y
  methods_on(dd, cnt, d)
}
## Exp 3 non-directional: arm B h = 1 + kappa (r - c)^2, c solved per slice so delta_exp = 0
one_rep_nd <- function(i, k, kappa) {
  d <- design_k(pool, k); d$Y <- NA_integer_; ddf <- numeric(0)
  for (s in levels(d$slice)) {
    j <- which(d$slice == s)
    if (d$C[j[1]] == "0") d$Y[j] <- sim_slice(d$r[j], 0, M_MED, SIZE)
    else { cc <- solve_c_nd(d$r[j], kappa); d$Y[j] <- sim_slice_nd(d$r[j], kappa, M_MED, SIZE)
           ddf <- c(ddf, nd_domain_diff(d$r[j], kappa, cc)) }
  }
  cnt <- sim_panel(d, c(0, 0), G_BG, M_MED, SIZE); cnt[1, ] <- d$Y
  out <- methods_on(d, cnt, d); out["I2"] <- mean(ddf); names(out)[names(out) == "I2"] <- "aux"; out
}
## Exp 6: G_FDR genes on ONE substrate per replicate; genes 1..G1 gained, rest null (they are the backdrop)
one_rep_fdr <- function(i, k) {
  d <- design_k(pool, k); TH1 <- theta_for("gained", 0.05)
  truth <- c(rep(TRUE, G1_FDR), rep(FALSE, G_FDR - G1_FDR))
  cnt <- matrix(0L, G_FDR, nrow(d), dimnames = list(paste0("g", seq_len(G_FDR)), paste0("c", seq_len(nrow(d)))))
  for (g in seq_len(G_FDR)) cnt[g, ] <- sim_design(d, if (truth[g]) TH1 else c(0, 0), M_MED, SIZE, tau_theta = 0.15)$Y
  pv <- mclapply(seq_len(G_FDR), function(g) { d$Y <- cnt[g, ]
    c(proposed = m_proposed(d)$pval, gam_fs = m_gam(d, "fs"), pbulk_de = m_pseudobulk_de(d)) }, mc.cores = CORES, mc.set.seed = TRUE)
  P <- do.call(rbind, pv); ds <- NULL; invisible(capture.output(ds <- unname(m_despace2_panel(cnt, d)[rownames(cnt)])))
  cbind(P, despace2 = ds)
}

run_cell <- function(id) {
  cl <- CELLS[id, ]
  key <- sprintf("exp%s|k%d|%s|D%.3f|tau%.2f", cl$exp, cl$k, cl$scen, cl$D, cl$tau)
  f <- file.path(OUTDIR, sprintf("cell_%02d.rds", id))
  if (file.exists(f)) { cat("skip", key, "\n"); return(invisible(NULL)) }
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260918 + id)
  cat("START", key, " cores =", CORES, "\n"); flush.console(); t0 <- Sys.time()
  if (cl$exp == "6") {
    reps <- vector("list", R_FDR)
    for (i in seq_len(R_FDR)) { reps[[i]] <- one_rep_fdr(i, cl$k); cat(sprintf("  fdr rep %d/%d (%.0f s)\n", i, R_FDR, as.numeric(Sys.time() - t0, units = "secs"))); flush.console() }
    out <- list(cell = cl, key = key, P = reps, truth = c(rep(TRUE, G1_FDR), rep(FALSE, G_FDR - G1_FDR)),
                G = G_FDR, G1 = G1_FDR, D = 0.05, tau = 0.15, M = M_MED, size = SIZE, n_rep = R_FDR)
  } else {
    res <- if (cl$scen %in% c("nd_moderate", "nd_strong"))
      mclapply(seq_len(R_N), one_rep_nd, k = cl$k, kappa = cl$kappa, mc.cores = CORES, mc.set.seed = TRUE)
    else mclapply(seq_len(R_N), one_rep_std, k = cl$k, theta_c = theta_for(cl$scen, cl$D),
                  Mv = c(M_MED, M_MED * cl$Mratio), tau = cl$tau, mc.cores = CORES, mc.set.seed = TRUE)
    ok <- vapply(res, function(z) is.numeric(z) && length(z) == 8, logical(1))
    P <- do.call(rbind, res[ok])
    out <- list(cell = cl, key = key, P = P, M = c(M_MED, M_MED * cl$Mratio), size = SIZE, n_rep = sum(ok))
  }
  out$secs <- as.numeric(Sys.time() - t0, units = "secs"); saveRDS(out, f)
  if (cl$exp != "6") { rj <- apply(out$P[, 1:5] < 0.05, 2, mean, na.rm = TRUE)
    cat(sprintf("DONE  %s  %d reps  %.0f s  reject@0.05: %s\n", key, out$n_rep, out$secs, paste(sprintf("%s=%.3f", names(rj), rj), collapse = " ")))
  } else cat(sprintf("DONE  %s  %d reps  %.0f s\n", key, out$n_rep, out$secs))
  flush.console()
}
args <- commandArgs(trailingOnly = TRUE)
for (id in if (length(args)) as.integer(args[1]) else CELLS$id) run_cell(id)
