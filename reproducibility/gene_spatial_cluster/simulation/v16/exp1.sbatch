#!/bin/bash
#SBATCH -J exp1_calib
#SBATCH -p shared
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 08:00:00
#SBATCH --array=1-12          # 12 cells: exp1 (k=10) null1/null2 x 4 tau = 8, exp1b (k=20) null1 x 4 tau = 4
#SBATCH -o logs/exp1_%a.out
#SBATCH -e logs/exp1_%a.err
# Plan v16 §0.4-0.5: Experiment 1 + 1b run FIRST and ALONE. Nothing else is submitted until
# summarise_exp1.R has been inspected (comparator Null 1 / Null 2 behaviour, c_{m,k,tau}, F0).
#
# Layout on FASRC (edit SUBSTRATE_ROOT if different):
#   $SUBSTRATE_ROOT/gene_spatial_cluster/gene_spatial_lineage_slices.csv
#   $SUBSTRATE_ROOT/output_threelayer/per_slice_morphogam/*_MU_xyt_mapped.csv
#   $SUBSTRATE_ROOT/gene_spatial_cluster/simulation/{pilot,v16}/   (this script runs from v16/)
# package layout: <root>/gene_spatial_cluster/simulation/v16 (here) -> root is three levels up
export SUBSTRATE_ROOT=${SUBSTRATE_ROOT:-$(cd "${SLURM_SUBMIT_DIR:-$PWD}/../../.." && pwd)}
export CORES=$SLURM_CPUS_PER_TASK
export R_N=${R_N:-200}
export G0=${G0:-50}
export G_GAM=${G_GAM:-50}       # GAM-fs genes scored per replicate; set 10 if the 8 h limit binds
mkdir -p logs exp1_out
source ../../config.sh 2>/dev/null || true
module purge 2>/dev/null; source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate rsim
# env check: the exact package set LOCKED_METHODS.md names
Rscript -e 'pk <- c("dplyr","readr","tidyr","purrr","metafor","mgcv","DESpace","SingleCellExperiment","edgeR");
  miss <- pk[!sapply(pk, requireNamespace, quietly=TRUE)];
  if (length(miss)) { cat("ENV FAIL, missing:", miss, "\n"); quit(status=1) }
  cat("R", R.version.string, "| DESpace", as.character(packageVersion("DESpace")), "| mgcv", as.character(packageVersion("mgcv")), "| metafor", as.character(packageVersion("metafor")), "\n")' || exit 1
test -f "$SUBSTRATE_ROOT/gene_spatial_cluster/gene_spatial_lineage_slices.csv" || { echo "SUBSTRATE_ROOT wrong: $SUBSTRATE_ROOT"; exit 1; }
Rscript run_exp1.R ${SLURM_ARRAY_TASK_ID}
