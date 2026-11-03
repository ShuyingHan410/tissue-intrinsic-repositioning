# Running the benchmark on FASRC — two designs × four methods

> This file documents the development-time harness (two within-Day 0 designs). The paper reports the
> Day 0 → Day 9 design, `../day0_day9/README_day0day9.md`.


Two benchmark designs, each scored the same way (self-match rank-1, pooled):
- **A "halfsplit"** — one Day0 slice → two stratified halves (10 slices).
- **B "slicepair"** — two DIFFERENT full Day0 slices → all C(10,2)=45 pairs.

Methods: your **three-step** (R, runs locally — cheap) + **moscot / PASTE2 / STalign**
(Python, on the cluster). All share `script/common/bench_common.py`, so each
design's input-cell filter (A = 7/11 types, B = 9/12 types) and pooling are applied
identically to every method.

## Layout
```
script/
  common/            bench_common.py, bench_heatmap.R, bench_proxy.R,
                     threestep_engine.R, run_bench.py, plot_results.R
  designA_halfsplit/ threestep.Rmd + moscot/paste2/stalign .Rmd (local A)
  designB_slicepair/ threestep.Rmd (local B; OT methods via run_bench.py --design slicepair)
  cluster/           setup_envs.sh, run_paste2.sbatch, this file
```

## 0. Paths / proxy
- Login `login.rc.fas.harvard.edu` (user `shuyhan`). `BENCH` = uploaded benchmark folder.
- **`adata.h5ad` is already on the cluster** at `/n/netscratch/rongma_lab/Lab/shuyhan/samples/adata.h5ad`
  — the jobs read it via `ADATA` (default set in the sbatch files). No need to re-upload it.
- Compute nodes have no internet → the jobs export `http://rcproxy.rc.fas.harvard.edu:3128`.

## 1. Upload the folder (from your Mac) — WITHOUT the 18 GB h5ad
```bash
RC_USER=shuyhan; RC_HOST=login.rc.fas.harvard.edu
RC_DEST=/n/home02/shuyhan/Rong_lab/benchmark
LOCAL="/Users/hanshuying/Desktop/PHchan/Lab/Ma lab/benchmark"
ssh ${RC_USER}@${RC_HOST} "mkdir -p ${RC_DEST}"
rsync -avhP --partial --inplace \
  --exclude '.DS_Store' --exclude '__pycache__' --exclude '.git' \
  --exclude 'adata.h5ad' \
  "${LOCAL}/" ${RC_USER}@${RC_HOST}:${RC_DEST}/
```
(adata.h5ad is skipped — it already lives at `/n/netscratch/rongma_lab/Lab/shuyhan/samples/`.
If your copy is elsewhere, set `ADATA=/path/to/adata.h5ad` when submitting.)

## 2. Build envs (once, interactive node)
Builds `otenv` (moscot+PASTE2) and `stalignenv` (STalign).
```bash
salloc -p test -c 4 --mem 16G -t 0-02:00
bash $RC_DEST/script/cluster/setup_envs.sh ; exit
```

## 3. Submit ALL OT methods × both designs
```bash
cd $RC_DEST/script/cluster
export BENCH=$RC_DEST
bash submit_all.sh            # moscot/paste2/stalign (CPU), A & B
squeue -u $USER
```
`submit_all.sh` loops `run_cpu.sbatch` (picks env by method: otenv for moscot/paste2,
stalignenv for stalign). Each writes
`benchmark_out/<design>/<method>/<method>_typematrix_{Fibroblast,Epithelial}.csv`.
Design B = 45 pairs × 2 lineages, so it runs longer (wall-time set in the sbatch).

Single combo instead of all:
```bash
METHOD=moscot DESIGN=slicepair sbatch run_cpu.sbatch
```

Each method also has a renderable Rmd (`designA_halfsplit/<method>.Rmd`,
`designB_slicepair/<method>.Rmd`) that calls the same `run_bench.run_method()` and
plots the heatmap into an HTML — handy for moscot/STalign locally
(`RETICULATE_PYTHON=<env>/bin/python Rscript -e 'rmarkdown::render("moscot.Rmd", knit_root_dir=getwd())'`).

## 4. Three-step method (your method) — run locally (fast, R)
```bash
cd "$LOCAL/script/designA_halfsplit" && Rscript -e 'rmarkdown::render("threestep.Rmd")'
cd "$LOCAL/script/designB_slicepair" && Rscript -e 'rmarkdown::render("threestep.Rmd")'
```
(No cluster needed — the KDE/JSD engine is cheap. It reads output_threelayer/.)

## 5. Pull results + plot everything (on your Mac)
```bash
rsync -avz ${RC_USER}@${RC_HOST}:${RC_DEST}/benchmark_out/ "${LOCAL}/benchmark_out/"
cd "${LOCAL}/script/common" && Rscript plot_results.R
# -> per-method heatmaps + benchmark_out/benchmark_summary.csv (designs × methods × lineages)
```

## Notes
- Raise/lower `--max_cells` in the sbatch files. Design B loads two full slices/pair.
- Verify partition names with `sinfo` (shared / gpu / bigmem).
