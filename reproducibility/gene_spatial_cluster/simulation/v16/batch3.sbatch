#!/bin/bash
#SBATCH -J batch3_robust
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 10:00:00
#SBATCH --array=1-85          # LIST=1 Rscript run_batch3.R prints the cell table
#SBATCH -o logs/batch3_%a.out
#SBATCH -e logs/batch3_%a.err
# Plan v16 §0.4 robustness experiments: 5 (real geometry + bootstrap), 7 (abundance), 7b (per-arm abundance),
# 8 (k = 20), 9 (deployed pipeline with lambda), 10 (Fibro 6), 11 (awkward shapes), 6b (FDR, correlated blocks).
# Submit after BATCH2_RESULTS.md is frozen. Most cells 5-25 min; Exp 9 cells ~30-60 min (300 blanks/rep);
# cell 85 (6b) is FDR-sized (~2 h).
export SUBSTRATE_ROOT=${SUBSTRATE_ROOT:-$(cd "${SLURM_SUBMIT_DIR:-$PWD}/../../.." && pwd)}
export CORES=$SLURM_CPUS_PER_TASK
export R_N=${R_N:-200}
export R_FDR=${R_FDR:-50}
mkdir -p logs batch3_out
module purge 2>/dev/null; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate rsim
Rscript -e 'pk <- c("dplyr","readr","tidyr","purrr","metafor","mgcv","DESpace","SingleCellExperiment","edgeR");
  miss <- pk[!sapply(pk, requireNamespace, quietly=TRUE)];
  if (length(miss)) { cat("ENV FAIL, missing:", miss, "\n"); quit(status=1) }
  cat("R", R.version.string, "| DESpace", as.character(packageVersion("DESpace")), "| mgcv", as.character(packageVersion("mgcv")), "| metafor", as.character(packageVersion("metafor")), "\n")' || exit 1
test -f "$SUBSTRATE_ROOT/gene_spatial_cluster/gene_spatial_lineage_slices.csv" || { echo "SUBSTRATE_ROOT wrong: $SUBSTRATE_ROOT"; exit 1; }
test -f ../pilot/null3_resampler.R || { echo "missing ../pilot/null3_resampler.R (copy from local simulation/pilot/)"; exit 1; }
test -f ../pilot/null3A_smoke.rds  || { echo "missing ../pilot/null3A_smoke.rds  (copy from local simulation/pilot/)"; exit 1; }
Rscript run_batch3.R ${SLURM_ARRAY_TASK_ID}
