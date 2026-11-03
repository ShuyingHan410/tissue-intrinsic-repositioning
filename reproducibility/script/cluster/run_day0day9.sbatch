#!/bin/bash
#SBATCH -J d0d9_arr
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem 16G
#SBATCH -t 0-04:00            # ONE Day0xDay9 pair (Fibroblast only) is small; 4 h is headroom
#SBATCH --array=0-249         # 10 Day0 x 25 Day9 (adata DSS9 slices) = 250 pairs; out-of-range units auto-skip
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
set -e
# Day0->Day9 CROSS-CONDITION benchmark, Fibroblast only. One SLURM array task ==
# one (Day0 Healthy slice, Day9 DSS9 slice) pair; writes
#   benchmark_out/day0_day9/<method>/partials/<method>_Fibroblast_unitNNN.csv
# Pool + score afterwards (see README_day0day9.md). CPU methods here (moscot/paste2/
# stalign/slat/moscot_ub).
#
# Submit e.g.:  METHOD=moscot sbatch run_day0day9.sbatch
BENCH=${BENCH:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}   # the reproducibility folder
ADATA=${ADATA:-$BENCH/adata.h5ad}
METHOD=${METHOD:?set METHOD=moscot|moscot_ub|paste2|stalign|slat}
MAXC=${MAXC:-4000}; [ "$METHOD" = paste2 ] && MAXC=${MAXC:-2500}

export http_proxy=http://rcproxy.rc.fas.harvard.edu:3128
export https_proxy=$http_proxy
CB="$(conda info --base 2>/dev/null)"; [ -d "$CB" ] || CB="${CONDA_BASE:?set CONDA_BASE to your conda installation}"
source "$CB/etc/profile.d/conda.sh"
[ -f "$CB/etc/profile.d/mamba.sh" ] && source "$CB/etc/profile.d/mamba.sh"
command -v mamba >/dev/null 2>&1 || mamba() { conda "$@"; }
case $METHOD in
  moscot|moscot_ub|paste2) mamba activate otenv ;;
  slat)          mamba activate slatenv ;;
  stalign)       mamba activate stalignenv ;;
  cast|cast_stack) mamba activate castenv ;;
  *) echo "bad METHOD $METHOD"; exit 1 ;;
esac
python "$BENCH/script/common/run_bench.py" \
    --method "$METHOD" --design day0day9 \
    --adata "$ADATA" \
    --out   "$BENCH/benchmark_out/day0_day9/$METHOD" \
    --lineages Fibroblast \
    --sample_type_a Healthy --sample_type_b DSS9 \
    --max_cells "$MAXC" \
    --unit "$SLURM_ARRAY_TASK_ID"
echo "DONE unit $SLURM_ARRAY_TASK_ID $METHOD day0day9 $(date)"
