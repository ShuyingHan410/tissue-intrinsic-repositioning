#!/bin/bash
#SBATCH -J gene_sweep_stratified
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:45:00
#SBATCH --array=1-33      # 990 cols / CHUNK_SIZE=30; check with: ls chunks/*.csv.gz | wc -l
#SBATCH -o logs/sweep_strat_%a.out
#SBATCH -e logs/sweep_strat_%a.err
source config.sh
eval "$ENV_SETUP"        # see config.sh: conda activate, or module load
# --- ENV CHECK: fail immediately with a useful message, not 40 lines into the script ---
command -v Rscript >/dev/null && Rscript -e 'q(status=as.integer(any(!sapply(c("dplyr","readr","tidyr","purrr","metafor"),requireNamespace,quietly=TRUE))))' 2>/dev/null || {
  echo "ENV FAIL: R or its packages missing on $(hostname)."
  echo "  which Rscript : $(which Rscript 2>&1)"
  echo "  Fix ENV_SETUP in config.sh (module load R with the SAME version you installed packages under)."
  exit 1; }
Rscript 01c_stratified_chunk.R ${SLURM_ARRAY_TASK_ID}
