#!/bin/bash
# Run FIRST on the cluster:  bash 99_check_paths.sh
# Verifies inputs, coordinate orientation, and deps BEFORE burning a job.
# Every check reports OK or FAIL - nothing is allowed to fail silently.
source config.sh
ok=1
fail() { echo "  FAIL $*"; ok=0; }

echo "== tools =="
command -v python3 >/dev/null && echo "  OK   python3  $(python3 --version 2>&1)" || fail "python3 not on PATH"
command -v Rscript >/dev/null && echo "  OK   Rscript  $(Rscript --version 2>&1 | head -1)" || \
  fail "Rscript not on PATH  ->  module load R    (or: module spider R)"

echo "== paths =="
if [ -f "$ADATA" ]; then echo "  OK   ADATA      $ADATA ($(du -h "$ADATA" | cut -f1))"
else fail "ADATA      $ADATA not found"
     [ -f "$ADATA_ALT1" ] && echo "       try instead: $ADATA_ALT1"; fi
if [ -d "$COORD_DIR" ]; then
  n=$(ls "$COORD_DIR"/*_MU_xyt_mapped.csv 2>/dev/null | wc -l)
  echo "  OK   COORD_DIR  $COORD_DIR  ($n MU files; expect 31)"
  [ "$n" -eq 0 ] && fail "no *_MU_xyt_mapped.csv in COORD_DIR"
else fail "COORD_DIR  $COORD_DIR not found"
     [ -d "$COORD_ALT1" ] && echo "       try instead: $COORD_ALT1"; fi
echo "  OUT   $OUT"

echo "== h5ad contents (must match the object the pipeline was built on) =="
if command -v python3 >/dev/null && python3 -c "import h5py" 2>/dev/null; then
  python3 - <<'PY' || fail "could not read $ADATA"
import h5py, os
f = h5py.File(os.environ["ADATA"], "r")
X = f["raw/X"]; var = f["raw/var/_index"]
nb = sum(1 for g in var[:] if (g.decode() if isinstance(g, bytes) else g).startswith("Blank"))
print(f"  {'OK  ' if X.shape==(1348582,990) else 'FAIL'} raw/X shape {X.shape}  (expect (1348582, 990))")
print(f"  {'OK  ' if nb==48 else 'WARN'} {nb} Blank barcodes  (pipeline assumes 48)")
for k in ("obs/Slice_ID","obs/Tier2","obs/Cell_ID"):
    print(f"  {'OK  ' if k in f else 'FAIL'} {k}")
PY
else echo "  SKIP (h5py unavailable - see dependencies below)"; fi

echo "== orientation (coords MUST be harmonised: lumen=0 -> outer=1) =="
if command -v Rscript >/dev/null; then
  out=$(Rscript -e '
d <- Sys.getenv("COORD_DIR"); fs <- list.files(d,"_(MU|SM|ME)_xyt_mapped[.]csv$",full.names=TRUE)
if (!length(fs)) { cat("FAIL no coordinate files\n"); quit(status=0) }
bad <- 0; amb <- 0
for (f in fs) { x <- try(read.csv(f), silent=TRUE); if (inherits(x,"try-error")) next
  x <- x[is.finite(x$x)&is.finite(x$y)&is.finite(x$r_for_jsd),]; if (nrow(x)<50) next
  cr <- suppressWarnings(cor(x$r_for_jsd, sqrt((x$x-median(x$x))^2+(x$y-median(x$y))^2), method="spearman"))
  if (is.finite(cr)) { if (cr < -0.05) bad <- bad+1 else if (abs(cr) <= 0.05) amb <- amb+1 } }
cat(sprintf("%s %d files | %d still INVERTED (want 0) | %d ambiguous (expected 5)\n",
            ifelse(bad==0,"OK  ","FAIL"), length(fs), bad, amb))' 2>&1 | grep -vE "UTF-8|^\[1\]")
  echo "  $out"
  echo "$out" | grep -q "^FAIL" && fail "coordinates are NOT harmonised - upload the post-2026-08-10 copy"
else fail "cannot verify orientation (no Rscript) - this is the check that matters most"; fi

echo "== dependencies =="
python3 -c "import h5py,numpy,pandas; print('  OK   python: h5py',h5py.__version__,'numpy',numpy.__version__,'pandas',pandas.__version__)" 2>/dev/null \
  || fail "python deps  ->  mamba install -c conda-forge h5py numpy pandas"
if command -v Rscript >/dev/null; then
  Rscript -e 'm <- c("dplyr","readr","tidyr"); miss <- m[!sapply(m,requireNamespace,quietly=TRUE)]
    if (length(miss)) cat("MISSING:", paste(miss,collapse=" "), "\n") else cat("  OK   R: dplyr readr tidyr\n")' 2>&1 \
    | grep -vE "UTF-8|^\[1\]" | sed 's/^MISSING/  FAIL R packages MISSING/'
  Rscript -e 'm<-c("dplyr","readr","tidyr"); quit(status=as.integer(any(!sapply(m,requireNamespace,quietly=TRUE))))' 2>/dev/null \
    || fail "install with:  Rscript -e 'install.packages(c(\"dplyr\",\"readr\",\"tidyr\"), repos=\"https://cloud.r-project.org\")'"
fi

echo "== chunk arithmetic =="
python3 -c "
import math,os
c=int(os.environ['CHUNK_SIZE']); print(f'  990 cols / {c} = {math.ceil(990/c)} chunks; array must be 1-{math.ceil(990/c)}')" 2>/dev/null \
  || echo "  SKIP (python unavailable)"
grep -o 'array=1-[0-9]*' 01_sweep.sbatch | sed 's/^/  currently: --/'

echo
[ "$ok" -eq 1 ] && echo "== ALL CHECKS PASSED - safe to run: bash submit_all.sh" \
                || echo "== FIX THE FAILURES ABOVE BEFORE SUBMITTING"
exit $((1-ok))
