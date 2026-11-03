#!/bin/bash
#SBATCH -J batch2_power
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH --array=1-30          # 1-24 Exp 2 (gained/reversed x 3 |D| x 4 tau), 25-28 Exp 3, 29-30 Exp 6 FDR (k=10, k=20)
#SBATCH -o logs/batch2_%a.out
#SBATCH -e logs/batch2_%a.err
# Plan v16 §0.5 step 4: submit ONLY after summarise_exp1.R has been run and the comparator
# Null 1 / Null 2 tables inspected. Cells 1-28 take minutes to ~40 min each; the two FDR
# cells take several hours (500 genes x 50 replicates, GAM-fs dominated).
export SUBSTRATE_ROOT=${SUBSTRATE_ROOT:-$(cd "${SLURM_SUBMIT_DIR:-$PWD}/../../.." && pwd)}
export CORES=$SLURM_CPUS_PER_TASK
export R_N=${R_N:-200}
export R_FDR=${R_FDR:-50}
mkdir -p logs batch2_out
module purge 2>/dev/null; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate rsim
Rscript -e 'pk <- c("dplyr","readr","tidyr","purrr","metafor","mgcv","DESpace","SingleCellExperiment","edgeR");
  miss <- pk[!sapply(pk, requireNamespace, quietly=TRUE)];
  if (length(miss)) { cat("ENV FAIL, missing:", miss, "\n"); quit(status=1) }
  cat("R", R.version.string, "| DESpace", as.character(packageVersion("DESpace")), "| mgcv", as.character(packageVersion("mgcv")), "| metafor", as.character(packageVersion("metafor")), "\n")' || exit 1
test -f "$SUBSTRATE_ROOT/gene_spatial_cluster/gene_spatial_lineage_slices.csv" || { echo "SUBSTRATE_ROOT wrong: $SUBSTRATE_ROOT"; exit 1; }
test -f exp1_thresholds.csv || echo "WARNING: exp1_thresholds.csv not found -- run summarise_exp1.R before interpreting batch 2"
Rscript run_batch2.R ${SLURM_ARRAY_TASK_ID}
