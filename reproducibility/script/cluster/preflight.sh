#!/bin/bash
# Pre-flight readiness check for submit_all.sh. Run on the LOGIN node:
#   bash preflight.sh
# Prints [ OK ] / [FAIL] / [WARN] per item and a final verdict.
BENCH=${BENCH:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}   # the reproducibility folder
ADATA=${ADATA:-$BENCH/adata.h5ad}
CM="$BENCH/script/common"
# DGL needs the conda env's newer libstdc++ (Rocky8 system one is too old)
export LD_LIBRARY_PATH="$(conda info --base 2>/dev/null)/envs/castenv/lib:$LD_LIBRARY_PATH"
FAILED=""
ok()   { printf "  [ OK ] %s\n" "$1"; }
bad()  { printf "  [FAIL] %s\n" "$1"; FAILED=1; }
warn() { printf "  [WARN] %s\n" "$1"; }

echo "=== 1. Files & data ==="
[ -f "$CM/run_bench.py" ]   && ok "run_bench.py"    || bad "run_bench.py missing ($CM)"
[ -f "$CM/bench_common.py" ]&& ok "bench_common.py" || bad "bench_common.py missing"
[ -f "$ADATA" ] && ok "adata.h5ad ($(du -h "$ADATA" 2>/dev/null | cut -f1))" || bad "adata.h5ad missing at $ADATA"
for f in run_cpu.sbatch run_cast.sbatch submit_all.sh; do
  [ -f "$BENCH/script/cluster/$f" ] && ok "$f" || bad "$f missing"
done

echo "=== 2. Conda envs import ==="
envok(){ conda run -n "$1" python -c "$2" >/dev/null 2>&1; }
envok otenv     "import moscot,scanpy,paste2"                    && ok "otenv (moscot+paste2)"   || bad "otenv import failed"
envok stalignenv "from STalign import STalign; import torch,sklearn,h5py,pandas" && ok "stalignenv (STalign)" || bad "stalignenv import failed"
if envok castenv "import torch,dgl; from CAST import CAST_MARK"; then ok "castenv (CAST+dgl)"
elif conda env list | awk '{print $1}' | grep -qx castenv;      then warn "castenv exists but full import failed on login (likely no GPU driver) - verify on a GPU node"
else bad "castenv missing"; fi

echo "=== 3. Driver reads the real data ==="
if conda run -n stalignenv python - <<PY >/dev/null 2>&1
import sys; sys.path.insert(0,"$CM"); import bench_common as bc, run_bench
p=dict(adata_path="$ADATA",sample_type="Healthy",lineage="Fibroblast",design="slicepair")
assert len(bc.all_slice_ids(p))==10
assert len(bc.all_units(p))==45
assert len(bc.type_set(p))==12
assert set(run_bench.METHODS)=={"moscot","paste2","stalign","cast"}
PY
then ok "10 slices / 45 pairs / 12 types / 4 methods"
else bad "driver could not read data (run the check manually to see the error)"; fi

echo
if [ -z "$FAILED" ]; then echo "READY  ->  bash submit_all.sh"
else echo "NOT READY - fix the [FAIL] items above (WARN on castenv is OK if the GPU job runs)"; fi
