#!/bin/bash
#SBATCH -J gene_extract
#SBATCH -p shared                # adjust to your FASRC partition
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o logs/extract_%j.out
#SBATCH -e logs/extract_%j.err
source config.sh
eval "$ENV_SETUP"        # see config.sh: conda activate, or module load
# --- ENV CHECK: fail immediately with a useful message, not 40 lines into the script ---
python3 -c "import h5py, numpy, pandas" 2>/dev/null || {
  echo "ENV FAIL: python deps missing on $(hostname)."
  echo "  CONDA_BASE_DIR=$CONDA_BASE_DIR"
  echo "  which python3 : $(which python3 2>&1)"
  echo "  Fix ENV_SETUP in config.sh so the compute node sees the same env as your login shell."
  exit 1; }
python3 00_extract_all.py
