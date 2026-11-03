# ---- paths -----------------------------------------------------------------------------
# Defaults assume the reproducibility layout (this folder = reproducibility/gene_spatial_cluster);
# export ADATA / COORD_DIR before sourcing to point elsewhere (e.g. on a cluster).
export ADATA=${ADATA:-$PWD/../adata.h5ad}
export COORD_DIR=${COORD_DIR:-$PWD/../output_threelayer/per_slice_morphogam}
export OUT=$PWD                                # this folder (chunks/, results/, logs/)

# Fallbacks if the layout differs - first existing path wins (checked by 99_check_paths.sh)
export ADATA_ALT1=${ADATA_ALT1:-$PWD/../../adata.h5ad}
export COORD_ALT1=${COORD_ALT1:-$PWD/../../output_threelayer/per_slice_morphogam}

# ---- fixed pipeline parameters (validated locally) -----------------------------------
export MIN_CELLS=30      # type must have >= this many cells in a slice
export MIN_MOL=20        # gene must have >= this many molecules in the type x slice
export MIN_SLICES=3      # >= this many eligible SLICES per side (slice-pair design)
export MIN_RATE_PER_1K=8 # ...OR >= this many molecules per 1000 cells (scales with type size)
export MIN_COVERAGE=0.5  # >= this fraction of available slices must be eligible
export CHUNK_SIZE=30     # genes per array task  (990 cols = 942 genes + 48 blanks -> 33 chunks)

# ---- how compute nodes get python + R -------------------------------------------------
# NOTE: `source ~/.bashrc` does NOT work in batch jobs - most .bashrc files return early
# when non-interactive, so conda never activates. Source conda's profile script instead.
#
# Find your base with:  conda info --base      (then hardcode it below if the guesses fail)
for _cb in "$CONDA_BASE_DIR" "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge" \
           "$HOME/miniforge3" "/n/sw/Mambaforge-22.11.1-4"; do
  if [ -f "$_cb/etc/profile.d/conda.sh" ]; then export CONDA_BASE_DIR="$_cb"; break; fi
done

# Runs on every compute node before python/R are used. Edit to match YOUR setup:
#   python from conda, R from the module system (what you have now)
export ENV_SETUP='
  if [ -n "$CONDA_BASE_DIR" ]; then . "$CONDA_BASE_DIR/etc/profile.d/conda.sh"; conda activate base; fi
  module load R 2>/dev/null || true
'
# Alternatives:
#   everything from modules:   export ENV_SETUP="module load python; module load R"
#   everything from conda:     export ENV_SETUP=". $CONDA_BASE_DIR/etc/profile.d/conda.sh; conda activate base"
