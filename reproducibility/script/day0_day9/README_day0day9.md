# Day0 → Day9 Fibroblast benchmark — MorphoGAM vs moscot / PASTE2 / SLAT / STalign


Cross-condition comparison: does each method recover the paper-supported Day0→Day9
fibroblast polarizations? Unlike the existing within-Day0 benchmark, every unit is a
**(Day0 Healthy slice, Day9 DSS9 slice)** pair, and the matrix is directional
(**rows = Day0 fibro types, cols = Day9 fibro/IAF types**).

## Ground truth (Day0 origin → Day9 target)
| origin | target(s) | compartment |
|---|---|---|
| Fibro 2 | IAF 4 | crypt top (MU) |
| Fibro 6 | IAF 2 **or** IAF 3 (best of the two) | crypt base (MU) |
| Fibro 7 | Fibro 4 | submucosa (SM/ME) |

## Metric (identical for all methods)
For each Day0 origin row, the **row-wise rank of its cognate Day9 target among that pair's
Day9 candidates** (rank 1 = best), pooled over all testable Day0×Day9 pairs → median rank,
mean normalized rank (0 = best), Top-1/3/5 %. Scored by `script/common/recovery.R`.

## Per-pair type sets
No fixed type list: each side of each pair keeps only its own Tier2 types with **≥30 cells**
(`day0day9_min_cells`), so the candidate set varies per pair — matching the MorphoGAM pipeline.

## Inputs
- **MorphoGAM** reads the per-slice coordinate CSVs `output_threelayer/per_slice_morphogam/`
  (manual). CSV-based, runs locally (cheap).
- **moscot / PASTE2 / SLAT / STalign** read `adata.h5ad` (Sample_type `Healthy` = Day0,
  `DSS9` = Day9), align a Day0 slice to a Day9 slice → cell coupling → type matrix. Run on FASRC.

## What gets generated
```
benchmark_out/day0_day9/
  morphogam/
    partials/morphogam_Fibroblast_unitNNN.csv     # per-pair Day0×Day9 similarity matrix
    morphogam_recovery_Fibroblast.csv             # pooled recovery (this method)
    morphogam_scored.rds
  moscot/  paste2/  slat/  stalign/
    partials/<method>_Fibroblast_unitNNN.csv       # per-pair type matrix (from alignment)
    <method>_typematrix_Fibroblast.csv             # pooled average
  recovery_comparison_Fibroblast.csv               # <- the head-to-head table (all methods)
```
Unit index NNN = position in the Day0×Day9 grid (`for d0 in sorted(Day0) for d9 in sorted(Day9)`),
identical across all methods.

---

## How to run

### 1. MorphoGAM (local, ~a few min)
```bash
cd script/day0_day9
Rscript -e 'rmarkdown::render("morphogam_recovery.Rmd", knit_root_dir=getwd())'
```
Writes `benchmark_out/day0_day9/morphogam/{partials,morphogam_recovery_Fibroblast.csv}`.

### 2. The four alignment methods (FASRC)
Prereqs (same as the existing benchmark): folder uploaded to `$BENCH`, envs built
(`setup_envs.sh` → otenv / stalignenv / slatenv), `adata.h5ad` at `$ADATA`. See
`script/cluster/README_cluster.md` §1–2.

**Submit one array job per method** (one SLURM task per Day0×Day9 pair):
```bash
cd $BENCH/script/cluster
export BENCH=$BENCH
METHOD=moscot  sbatch run_day0day9.sbatch      # otenv  (CPU)
METHOD=paste2  MAXC=2500 sbatch run_day0day9.sbatch
METHOD=stalign sbatch run_day0day9.sbatch      # stalignenv (CPU)
METHOD=slat    sbatch run_day0day9.sbatch      # slatenv (CPU; see the SLAT section below)
squeue -u $USER
```
Each task writes one `partials/<method>_Fibroblast_unitNNN.csv`. The array is `0-209`
(≤10×21 pairs); out-of-range units auto-skip.

**Pool each method's partials** (after its array finishes — quick, no GPU):
```bash
for M in moscot paste2 slat stalign; do
  python $BENCH/script/common/run_bench.py --method $M --design day0day9 \
     --adata $ADATA --out $BENCH/benchmark_out/day0_day9/$M \
     --lineages Fibroblast --sample_type_a Healthy --sample_type_b DSS9 --pool
done
```

### 3. Pull results back (your Mac)
```bash
rsync -avz shuyhan@login.rc.fas.harvard.edu:$BENCH/benchmark_out/day0_day9/ \
      "/Users/hanshuying/Desktop/PHchan/Lab/Ma lab/benchmark/benchmark_out/day0_day9/"
```

### 4. Head-to-head comparison (local)
```bash
cd script/day0_day9
Rscript -e 'rmarkdown::render("recovery_compare.Rmd", knit_root_dir=getwd())'
```
Reads every method's `partials/`, scores them identically, and writes
`benchmark_out/day0_day9/recovery_comparison_Fibroblast.csv` + the ranked-method HTML.

---

## Notes
- Cross-condition alignment (Healthy→DSS9) is harder than Day0↔Day0; expect the alignment
  methods to be noisier than in the within-Day0 benchmark. MorphoGAM (distribution-based, no
  alignment) is not affected by that.
- `--max_cells` caps cells per side per pair; lower it if PASTE2/STalign are slow.
- The within-condition ordinal CI (`*_ci_summary`) is **skipped** for day0day9 — the
  rank-recovery in step 4 is the metric here.
- To use the improved coordinates instead of manual, point `per_slice_xyt_dir` in
  `morphogam_recovery.Rmd` at `per_slice_morphogam_NEW` (one line).

---
## SLAT and unbalanced moscot (added 2026-10-28)

Two comparators designed for pairs that differ (different conditions, partial overlap, changed
composition):

| METHOD | env | what it does |
|---|---|---|
| `slat` | slatenv | SLAT graph-adversarial alignment; soft top-20 embedding matching per Day9 cell (no shared-frame assumption) |
| `moscot_ub` | otenv | moscot fused GW with unbalanced marginals (tau 0.9): populations may appear / vanish |

```bash
sbatch setup_envs.sh                              # once: builds slatenv (see notes inside)
cd $BENCH/script/cluster; export BENCH=$BENCH
METHOD=slat      sbatch run_day0day9.sbatch
METHOD=moscot_ub sbatch run_day0day9.sbatch
```
Then pool and score exactly as for the other methods. Scoring note: comparator matrices are
coupling MASS and track target abundance; score them column-normalised (mass per unit of
target mass) as well as raw. `benchmark_out/day0_day9/recovery_comparison_raw_vs_abundance_corrected.csv`
holds both scorings for the existing methods; the expression-only baseline is in
`benchmark_out/day0_day9/expression/`.

---
## Full rerun procedure (2026-10-29): slat, moscot_ub

### 0. Upload (Mac). Only five files changed; the folder and adata.h5ad are already on FASRC
```bash
bash script/cluster/upload_day0day9.sh
```
(run_bench.py, bench_common.py, run_day0day9.sbatch, setup_envs.sh, this README.)

### 1. Build slatenv (once). LOGIN node, because PyG wheels come from data.pyg.org
```bash
ssh shuyhan@login.rc.fas.harvard.edu
source ~/miniconda3/etc/profile.d/conda.sh
mamba create -y -c conda-forge --override-channels -n slatenv python=3.11
mamba run -n slatenv pip install scSLAT scanpy anndata h5py faiss-cpu
mamba run -n slatenv install_pyg_dependencies
mamba run -n slatenv python -c "from scSLAT.model import Cal_Spatial_Net, load_anndatas, run_SLAT, spatial_match; print('slat ok')"
```
otenv (moscot_ub) already exists.

### 2. Smoke-test ONE unit per method before the arrays
```bash
cd /n/home02/shuyhan/Rong_lab/benchmark/script/cluster
METHOD=slat       sbatch --array=0 run_day0day9.sbatch
METHOD=moscot_ub  sbatch --array=0 run_day0day9.sbatch
squeue -u $USER
```
Pass = the .out ends with `DONE unit 0 ...` and
`benchmark_out/day0_day9/<method>/partials/<method>_Fibroblast_unit000.csv` exists with a
non-empty matrix (`head -3` it). Any Python traceback is in `d0d9_arr_<jobid>_0.err`.

### 3. Full arrays (250 units each; out-of-range units auto-skip)
```bash
METHOD=slat       sbatch run_day0day9.sbatch
METHOD=moscot_ub  sbatch run_day0day9.sbatch
```
Check completeness: `ls benchmark_out/day0_day9/<method>/partials | wc -l` -> 250.
Re-run missing units with `--array=<n1>,<n2>,...`.

### 4. Pull back (Mac) and rescore everything
```bash
rsync -avz shuyhan@login.rc.fas.harvard.edu:/n/home02/shuyhan/Rong_lab/benchmark/benchmark_out/day0_day9/ \
      "/Users/hanshuying/Desktop/PHchan/Lab/Ma lab/benchmark/benchmark_out/day0_day9/"
cd "/Users/hanshuying/Desktop/PHchan/Lab/Ma lab/benchmark"
Rscript script/day0_day9/rescore_all_methods.R
```
Writes `recovery_comparison_all_methods.csv` (per origin), `recovery_pooled_all_methods.csv`
(one row per method x scoring) and `per_pair_scores_all_methods.csv`. Report the
**column-normalised** scoring; use the ties-averaged variant for any hard-matching method.
