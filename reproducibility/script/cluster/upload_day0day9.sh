#!/bin/bash
# Upload ONLY the files changed for the Day0->Day9 rerun (slat / moscot_ub).
# The benchmark folder and adata.h5ad are already on FASRC; nothing else needs to move.
# Run from your Mac:   bash script/cluster/upload_day0day9.sh
set -e
RC=${RC:-user@login.cluster}            # set RC, DEST and LOCAL for your own machine and cluster
DEST=${DEST:-/path/on/cluster/benchmark}
LOCAL=${LOCAL:-$(pwd)}
cd "$LOCAL"
rsync -avhR \
  script/common/run_bench.py \
  script/common/bench_common.py \
  script/cluster/run_day0day9.sbatch \
  script/cluster/setup_envs.sh \
  script/day0_day9/README_day0day9.md \
  "$RC:$DEST/"
echo "uploaded to $RC:$DEST"
