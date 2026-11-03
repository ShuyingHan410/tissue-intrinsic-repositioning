## ============================================================
## Experiment 1 (k = 10: Null 1 + Null 2) and 1b (k = 20: Null 1)
## Plan v16 §0.4–0.5. Reference design: M = 64, NB size 1, matched pseudo-arms
## redrawn per replicate, gate SE0 < 0.05 & cells >= 30 (engine), R = 200.
## 50 null genes per replicate, scored for EVERY method (§0a A2). Saves raw
## p-value matrices (R x G0) per method so c_{m,k,tau} and F0_{m,k,tau} can be
## computed downstream (summarise_exp1.R).
##
## Usage (one cell per SLURM array task):
##   Rscript run_exp1.R <cell_id>         cell_id in 1..12, see CELLS below
##   Rscript run_exp1.R                   runs every cell sequentially (laptop)
## Env: R_N (200), G0 (50), G_GAM (50; GAM-fs genes scored per rep, <= G0),
##      CORES (SLURM_CPUS_PER_TASK or detectCores()-1), SUBSTRATE_ROOT
## ============================================================
suppressMessages({source("engine_v16.R"); library(DESpace); library(SingleCellExperiment); library(parallel)})

R_N   <- as.integer(Sys.getenv("R_N", "200"))
G0    <- as.integer(Sys.getenv("G0", "50"))
G_GAM <- min(G0, as.integer(Sys.getenv("G_GAM", "50")))
M_MED <- 64; SIZE <- 1
CORES <- as.integer(Sys.getenv("CORES", Sys.getenv("SLURM_CPUS_PER_TASK", max(1L, detectCores() - 1L))))
ROOT  <- Sys.getenv("SUBSTRATE_ROOT", "../../..")
OUTDIR <- Sys.getenv("EXP1_OUT", "exp1_out"); dir.create(OUTDIR, showWarnings = FALSE)

TAUS  <- c(0, 0.05, 0.15, 0.30)
CELLS <- rbind(expand.grid(exp = "1",  k = 10, null = c("null1", "null2"), tau = TAUS, stringsAsFactors = FALSE),
               expand.grid(exp = "1b", k = 20, null = "null1",             tau = TAUS, stringsAsFactors = FALSE))
CELLS$id <- seq_len(nrow(CELLS))

x <- load_substrate(ROOT); pool <- substrate_pool(x)
TH05 <- solve_theta(pool$r, 0.05)           # Null 2: positioned identically in both arms

one_rep <- function(i, k, theta_c, tau) {
  d   <- design_k(pool, k)
  cnt <- sim_panel(d, theta_c, G0, M_MED, SIZE, tau_theta = tau)      # G0 null genes, same theta both arms
  P <- matrix(NA_real_, G0, 5, dimnames = list(rownames(cnt), c("proposed", "gam_naive", "gam_fs", "despace2", "pbulk_de")))
  invisible(capture.output(P[, "despace2"] <- m_despace2_panel(cnt, d)[rownames(cnt)]))  # one call returns all G0; DESpace prints its design matrix
  for (g in seq_len(G0)) {
    d$Y <- cnt[g, ]
    P[g, "proposed"] <- m_proposed(d)$pval
    P[g, "pbulk_de"] <- m_pseudobulk_de(d)
    if (g <= G_GAM) { P[g, "gam_fs"] <- m_gam(d, "fs"); P[g, "gam_naive"] <- m_gam(d, "naive") }
  }
  P
}

run_cell <- function(id) {
  cl <- CELLS[id, ]
  key <- sprintf("exp%s|k%d|%s|tau%.2f", cl$exp, cl$k, cl$null, cl$tau)
  f <- file.path(OUTDIR, sprintf("cell_%02d.rds", id))
  if (file.exists(f)) { cat("skip", key, "\n"); return(invisible(NULL)) }
  theta_c <- if (cl$null == "null1") c(0, 0) else c(TH05, TH05)
  RNGkind("L'Ecuyer-CMRG"); set.seed(20260917 + id)
  cat("START", key, " R =", R_N, " G0 =", G0, " G_GAM =", G_GAM, " cores =", CORES, "\n"); flush.console()
  t0 <- Sys.time()
  res <- mclapply(seq_len(R_N), one_rep, k = cl$k, theta_c = theta_c, tau = cl$tau, mc.cores = CORES, mc.set.seed = TRUE)
  ok  <- vapply(res, is.matrix, logical(1))
  A   <- simplify2array(res[ok])                                   # G0 x 5 x R
  out <- list(cell = cl, key = key, theta_c = theta_c, M = M_MED, size = SIZE, G0 = G0, G_GAM = G_GAM,
              P = A, n_rep = sum(ok), secs = as.numeric(Sys.time() - t0, units = "secs"),
              session = sessionInfo()$otherPkgs |> names())
  saveRDS(out, f)
  t1 <- apply(A < 0.05, 2, mean, na.rm = TRUE)
  cat(sprintf("DONE  %s  %d reps  %.0f s  TypeI: %s\n", key, sum(ok), out$secs,
              paste(sprintf("%s=%.3f", names(t1), t1), collapse = " "))); flush.console()
}

args <- commandArgs(trailingOnly = TRUE)
ids <- if (length(args)) as.integer(args[1]) else CELLS$id
for (id in ids) run_cell(id)
