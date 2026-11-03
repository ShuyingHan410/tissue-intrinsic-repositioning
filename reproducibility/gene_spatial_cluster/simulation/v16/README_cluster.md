# v16 cluster stage — Experiment 1 + 1b first, alone

Plan: `../revised_gene_repositioning_simulation_plan_v16.md` §0 (operative). Local calibration stage:
`CALIBRATION_STAGE.md`. Gate frozen 2026-09-13 (`../../GATE_FREEZE_20260913.md`).

## Starting from an empty cluster account
Everything needed is in `cluster_package_v16.tar.gz` (43 MB; 187 MB unpacked), built 2026-09-13 from the
harmonised coordinates (inner = 0 / outer = 1; 26/27 evaluable slices verified, 0 inverted — see `MANIFEST.md`).
```
scp cluster_package_v16.tar.gz <user>@login.rc.fas.harvard.edu:/n/home02/shuyhan/Rong_lab/
ssh ... ; cd /n/home02/shuyhan/Rong_lab && tar -xzf cluster_package_v16.tar.gz
cd cluster_package_v16/gene_spatial_cluster/simulation/v16
sbatch exp1.sbatch
```
`exp1.sbatch` finds the substrate three levels above itself, so the tree must stay intact. Nothing else
(the real pipeline, chunks, adata) is required for the simulation.

### Layout inside the package
```
gene_spatial_cluster/simulation/pilot/          engine.R + LOCKED_METHODS.md (frozen comparators)
gene_spatial_cluster/simulation/v16/            this folder
gene_spatial_cluster/gene_spatial_lineage_slices.csv
output_threelayer/per_slice_morphogam/*_MU_xyt_mapped.csv
```
`engine.R::load_substrate(root)` expects `root/gene_spatial_cluster/...` and `root/output_threelayer/...`.
`exp1.sbatch` sets `SUBSTRATE_ROOT` to the package root automatically; override with
`export SUBSTRATE_ROOT=...` only if you move the pieces apart.

## R packages needed on the cluster (LOCKED_METHODS.md environment)
dplyr readr tidyr purrr metafor mgcv DESpace (2.4.x) SingleCellExperiment edgeR ggplot2
`exp1.sbatch` checks these and exits with a message if any is missing.

## Submit
```
cd .../gene_spatial_cluster/simulation/v16
sbatch exp1.sbatch            # array 1-12; ~2 h (k=10) to ~4 h (k=20) per cell at 8 cores
```
If the 8 h limit binds, resubmit with `G_GAM=10 sbatch --export=ALL,G_GAM=10 exp1.sbatch`
(GAM-fs scores 10 of the 50 null genes per replicate; proposed and DESpace2 still score all 50) —
and state that in the results.

## Inspect before anything else is submitted
```
Rscript summarise_exp1.R      # needs all 12 exp1_out/cell_*.rds
```
Outputs: `exp1_typeI.csv` (Type I per method x design, pass flag [3%,7%]), `exp1_thresholds.csv`
(design-specific c_{m,k,tau}), `exp1_F0.rds` (empirical null CDFs for the FDR repair), `exp1_null2.csv`
(Null 2 rejection rates — the "is there a gradient?" trap), `exp1_qq.png`.

Decision rule (§0.5): a comparator inside [3%, 7%] at a design uses nominal p < 0.05 there; outside, it
uses c_{m,k,tau} and its FDR uses BH on F0-calibrated p. A high Null 2 rate is a specification problem
(ordered-factor trap), not a threshold problem — diagnose before running power. Only then: Experiments
2, 3, 4B, 6; freeze; then 5, 7, 7b, 8, 9, 10, 11, 6b.

## Batch 3 — robustness (Experiments 5, 7, 7b, 8, 9, 10, 11, 6b); submit after BATCH2_RESULTS.md is frozen
Files to copy into the cluster `v16/`: `run_batch3.R`, `batch3.sbatch`, `summarise_batch3.R`.
Also needed in the cluster `../pilot/`: `null3_resampler.R` and `null3A_smoke.rds` (the frozen Null 3 tilt) — copy from
local `simulation/pilot/` if absent.
```
LIST=1 Rscript run_batch3.R   # prints the 85-cell table
sbatch batch3.sbatch          # array 1-85; most cells 5-25 min, Exp 9 cells 30-60 min, cell 85 (6b) ~2 h
Rscript summarise_batch3.R    # after all 85 exist; needs exp1_thresholds.csv, exp1_F0.rds, batch2_power.csv
```
Thresholds: pseudo-arm cells at M = 64 use Experiment 1's c_{m,k,tau}; every other design (real geometry, bootstrap,
M = 27/331, per-arm M, Fibro 6, mouse effects) uses its own in-batch Null 1 partner cell, so nominal and matched are
both reported everywhere. Exp 9 reports proposed with and without lambda plus the blank diagnostics.
Exp 10 (Fibro 6) runs at NB size 2 and M = 43 (64 x 28/42, the same all-panel measurement on both substrates).

## Local pieces already done (do not rerun)
Exp 0 (Fibro 2: size 1; Fibro 6: size 2), proposed-method Null 1 at k in {5,10,20} and tau in {0,.05,.15,.30},
tau -> I2 map, sigma_mouse estimate, gate dry run + freeze, engine smoke tests, proposed-only power previews.
