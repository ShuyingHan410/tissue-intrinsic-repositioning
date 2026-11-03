#!/bin/bash
#SBATCH -J gene_merge
#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:20:00
#SBATCH -o logs/merge_%j.out
#SBATCH -e logs/merge_%j.err
source config.sh
eval "$ENV_SETUP"        # see config.sh: conda activate, or module load
# --- ENV CHECK: fail immediately with a useful message, not 40 lines into the script ---
command -v Rscript >/dev/null && Rscript -e 'q(status=as.integer(any(!sapply(c("dplyr","readr","tidyr"),requireNamespace,quietly=TRUE))))' 2>/dev/null || {
  echo "ENV FAIL: R or its packages missing on $(hostname)."
  echo "  which Rscript : $(which Rscript 2>&1)"
  echo "  Fix ENV_SETUP in config.sh (module load R with the SAME version you installed packages under)."
  exit 1; }
Rscript 02_merge.R
