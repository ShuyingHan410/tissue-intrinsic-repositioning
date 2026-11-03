# Reproducibility: analyses of the paper

Code only. This folder holds every script behind the analyses in

> Shuying Han and Rong Ma. *Quantifying Cell and Gene Repositioning Across Tissue States Using
> Tissue-Intrinsic Geometric Coordinates.* Submitted to RECOMB 2027.

The method itself is the R package in the repository root (`../README.md`). This folder keeps the
directory layout of the original analysis so that every script's relative paths work unchanged;
**run all commands from this folder** unless a line says otherwise. Intermediate outputs are not in
git; they come either from the checkpoint release (Level 2) or from re-running the stages (Level 3).

There are three levels of reproduction, from minutes to weeks:

| Level | What you get | Needs | Time |
|---|---|---|---|
| **1. Method demonstration** | the package, its unit tests and the vignette on shipped example data | R | minutes |
| **2. Manuscript figures** | Figures 1–5 and the supplementary figures regenerated from documented checkpoints | R, Python, the checkpoint download (≈ 300 MB) | ≈ 15 min |
| **3. Full analysis from raw data** | every checkpoint rebuilt from the public atlas | the atlas (≈ 6 GB), MorphoGAM, a SLURM cluster for the comparators and the gene sweep; one curated step cannot be automated | days to weeks |

`../FIGURE_REPRODUCTION_MANIFEST.md` maps every submitted figure to its script, inputs, outputs and
verification status. `../RECOMB_REPRODUCIBILITY_AUDIT.md` lists which checks were run and passed.

---

## Level 1 — method demonstration (no download)

```r
remotes::install_github("ShuyingHan410/tissue-intrinsic-repositioning", build_vignettes = TRUE)
library(TICR); vignette("quick_start", package = "TICR")        # both analyses on shipped example data
testthat::test_dir(system.file("tests", package = "TICR"))        # or, from a clone: devtools::test()
```

The vignette runs the population correspondence on one Day 0 × Day 9 colon section pair and the
*Igfbp5* displacement of Figure 5D on all sections of Fibro 2 and IAF 4 (`inst/extdata/`).

---

## Level 2 — regenerate the manuscript figures from checkpoints

The checkpoint release contains the 93 curated coordinate files (the one input that cannot be rebuilt
automatically, see Stage 1) and the small result tables of stages 2–4. **Download location:** to be
filled in once the release is published (GitHub release of this repository or Zenodo; see the release
README for the file list and SHA-256 sums).

1. Unpack into this folder so that these paths exist:

```
reproducibility/
  output_threelayer/per_slice_morphogam/*.csv        93 coordinate files            (coordinates/)
  atlas_sections.csv                                 section table                  (checkpoints/figure_inputs/)
  benchmark_out/day0_day9/per_pair_scores_all_methods.csv, recovery_*.csv          (checkpoints/benchmark/)
  benchmark_out/day0_day9/morphogam/partials/*.csv   210 per-pair matrices          (checkpoints/benchmark/morphogam_partials/)
  benchmark_out/day9_type_counts_full.csv, day9_fibro_totals.csv                   (checkpoints/figure_inputs/benchmark_out/)
  gene_spatial_cluster/contrastsC.csv, armsC.csv, lambdaC.csv, blankC_by_layer.csv (checkpoints/gene_displacement/)
  gene_spatial_cluster/simulation/v16/*.csv, v17_revision/*.csv                     (checkpoints/simulation/)
  manuscript_figures/fig1_assets/*, fig23_assets/*                                  (checkpoints/figure_inputs/)
```

2. Regenerate (each script writes its figure next to itself; Python needs `numpy pandas h5py matplotlib`,
R needs `ggplot2 patchwork dplyr readr tidyr`):

```bash
cd manuscript_figures
python3 make_fig2_correspondence.py            # Fig. 2  (also writes fig23_assets/fig2A_pooled_similarity.csv)
python3 make_fig1_v11.py "$PWD"                # Fig. 1
python3 make_fig3_benchmark_v2.py              # Fig. 3
Rscript  make_fig4_simulation.R                # Fig. 4
python3 make_fig5_repositioning_v2.py          # Fig. 5
Rscript  make_fig3_gene_simulation.R           # supplement: FDR and ablation panels (figS_fdr_ablation)
python3 make_figS_benchmark_sensitivity.py     # supplement: benchmark scoring sensitivity
Rscript  make_figS_skin_validation.R           # supplement: skin application (values from the skin report)
```

The figure-input tables under `fig1_assets/` and `fig23_assets/` and the two `day9_*` tables are
themselves derived; `make_figure_assets.py` rebuilds all of them from the stage outputs plus the atlas
(Level 3) and was verified to reproduce the shipped files exactly. Stage-3 tables can be re-derived
from `contrastsC.csv` with `gene_spatial_cluster/04_figure5_calls.R`.

PDF outputs differ from the submitted PDFs only in embedded creation timestamps; the PNG outputs
are byte-identical to the submitted ones (see the manifest).

---

## Level 3 — full analysis from the public atlas

Legend: **[auto]** runs unattended; **[cluster]** SLURM array jobs, hours to days; **[manual]** needs
the user's intervention; **[checkpoint]** output that is distributed in Level 2.

### Data

| Dataset | Source | Used for |
|---|---|---|
| Mouse DSS-colitis colon MERFISH (Cadinu et al., *Cell* 2024) | Dryad [doi:10.5061/dryad.rjdfn2zh3](https://doi.org/10.5061/dryad.rjdfn2zh3), CC0 (`adata.h5ad`, 940 genes + 50 blank barcodes, 1.35 M cells, 52 sections); also Bioconductor `MerfishData::MouseColonIbdCadinu2024()` | 10 Day 0 and 21 Day 9 sections |
| Mouse MC903 / oxazolone skin MERFISH (Imanishi et al., *Cell Reports* 2025) | Mendeley Data [10.17632/s3hgv4fp7c.1](https://doi.org/10.17632/s3hgv4fp7c.1) (Seurat v5 object) | 16 full-thickness strips |

Three small files that the pipelines read and that are not derivable from code are in git:
`gene_spatial_cluster/genes_index.csv` (gene → chunk map with blank flags),
`gene_spatial_cluster/tableS2_expressed.csv` (the atlas's Table S2 expression calls) and
`gene_spatial_cluster/simulation/pilot/null3A_smoke.rds` (the frozen moment-balanced tilt that
defines the whole-population-movement scenario, produced with `pilot/null3_resampler.R`).

### Step 1 — section tables **[auto]**

```bash
# adata.h5ad downloaded from Dryad into this folder
python3 export_colon_sections.py adata.h5ad all_slices     # -> all_slices/<Slice_ID>.csv (52 files) and atlas_sections.csv
```

One CSV per section with the atlas `obs` columns (`x, y, Mouse_ID, Technical_repeat_number, Sample_type,
Slice_ID, FOV, Cell_ID, sample, N_genes, Tier1, Tier2, Tier3, Leiden_neigh` and the UMAP columns);
verified identical to the tables used in the paper for all 52 sections.

### Software environment

R ≥ 4.5 (analyses were run with R 4.6.1 locally and R 4.5.3 on the cluster):

```r
source("environment/install_R_packages.R")   # CRAN packages, MorphoGAM at the commit used, DESpace/edgeR for simulation comparators
```

Python ≥ 3.9 with `numpy pandas h5py matplotlib` (`environment/requirements.txt`) for the h5ad
extraction and the figure scripts. The cell-matching comparators need their own conda
environments (moscot + PASTE2, STalign, scSLAT); `script/cluster/setup_envs.sh`
builds them on a SLURM cluster and `script/cluster/README_cluster.md` documents the runs.

---


### Stage 1 — tissue-intrinsic coordinates

### Colon

**User intervention required.** 74 of the 93 section × layer fits of the paper are hand-traced
(`CurveFinderInteractive`), 13 are automatic and 6 automatic with `knn = 50` (see the manifest shipped
with the coordinates). Running this notebook therefore reproduces the automatic fits but not the
hand-traced ones; the coordinates used in the paper are distributed as a checkpoint (Level 2).


`three_layer.Rmd` produces `output_threelayer/per_slice_morphogam/`. For every section and layer
(`MU`, `SM`, `ME` from the atlas neighbourhood labels; populations with < 30 cells in a layer are
dropped) it calls MorphoGAM `CurveFinder` and, where the default kNN graph is disconnected, retries
with `knn = 50` and finally `CurveFinderInteractive` (hand-traced medial curve; the `morphogam_manual`
and `morphogam_knn` columns record which rung each file came from). The notebook reads the section tables from `all_slices/`
(or `$ALL_SLICES`), written by `export_colon_sections.py` (Level 3, step 1).

`Rscript output_threelayer/flip_orientation.R` then harmonises orientation so that `r_for_jsd` runs
from the lumen (0) to the outer boundary (1) in every file, using the Spearman correlation between
`r` and distance from the layer centroid, and writes an audit table.

`morpho_improve/final/curve_finder.R` holds the MorphoGAM helpers (connectivity guard, projection
of all cells onto a curve fitted on a subsample, orientation rule) that the skin fitter sources; the
package function `fit_layer_coordinates()` wraps the same code.

### Skin

```bash
cd Skin_MC903_MERFISH
Rscript code/extract_cells.R          # Seurat object -> data/mc903_cells.csv.gz
Rscript code/export_counts.R          # 500-gene counts -> data/mc903_counts.csv.gz
Rscript code/export_slices.R          # depth compartments (code/skin_layers.R) -> all_slices/<strip>.csv
Rscript code/run_plain_morphogam.R    # one CurveFinder call per strip x compartment -> output/morphogam_auto/
Rscript code/harmonise_orientation.R  # surface 0 -> deep 1
```

Compartments come from anatomical landmarks only (`skin_layers.R`): distance to the connected
keratinocyte sheet, the first adipocytes and the panniculus; the papillary dermis is the upper 40 %
of the dermis.

---


### Stage 2 — cell-population correspondence and benchmark

`script/common/threestep_engine.R` implements the similarity score of the paper: Gaussian-KDE
densities of `r_for_jsd` on [0, 1] per layer (`bw = 2 × bw.nrd0`, floor 0.05, 512-point grid),
Jensen–Shannon distance, `S = exp(−D/σ)` with `σ` the median distance in the layer, layer weights
`sqrt(N_a N_b)` and the occupancy penalty `O_ab = Σ_h sqrt(π_a^h π_b^h)`. Populations need ≥ 30 cells
per evaluated layer.

```bash
Rscript -e 'rmarkdown::render("script/spatial_similarity_full_pipeline.Rmd")'   # method walk-through + application
Rscript -e 'rmarkdown::render("script/day0_day9/morphogam_recovery.Rmd")'       # all 10 x 21 pairs -> benchmark_out/day0_day9/morphogam/
```

Comparators run on the same section pairs through `script/common/run_bench.py --design day0day9
--method {moscot,paste2,slat,stalign}` (one SLURM array task per pair, `script/cluster/run_day0day9.sbatch`),
each aggregating its cell-level coupling to a population × population matrix under
`benchmark_out/day0_day9/<method>/partials/`. All methods are scored identically by
`script/common/recovery.R` (row rank of the atlas-supported target among the pair's Day 9 candidates;
ties take the best rank):

```bash
Rscript script/day0_day9/rescore_all_methods.R      # -> benchmark_out/day0_day9/recovery_*.csv, per_pair_scores_all_methods.csv
```

---


### Stage 3 — population-relative gene displacement (colitis atlas)

The analysis reported in the paper is stratified by anatomical layer: every correspondence × layer is
a self-contained analysis (δ per gene × population × layer × section, REML + Knapp–Hartung pooling
over the sections of that layer, Welch contrast Δ = δ(Day 9) − δ(Day 0)), and Benjamini–Hochberg runs
within each correspondence × layer family.

```bash
cd gene_spatial_cluster && source config.sh     # ADATA, COORD_DIR, thresholds
bash 99_check_paths.sh                          # verifies inputs and that every coordinate file is lumen = 0
python3 00_extract_all.py                       # h5ad -> meta.csv.gz (+ blank counts) and 33 gene chunks
python3 00b_blank_thresholds.py                 # per-barcode blank thresholds -> blank_thresholds.csv
sbatch 01c_stratified.sbatch                    # = Rscript 01c_stratified_chunk.R <chunk> for chunks 1-33:
                                                #   per gene x correspondence x layer x section delta and SE0, eligibility
                                                #   gate, REML + Knapp-Hartung per arm within each layer -> resultsC/part_<chunk>.csv
Rscript 02d_stratified_merge.R                  # bind, routing (D1_full / daughter_only / parent_only), Welch contrast Delta,
                                                #   blank-lambda diagnostics (lambdaC.csv) -> contrastsC.csv (2,030 rows; 892 D1_full)
Rscript 04_figure5_calls.R                      # BH on the raw Welch p within each of the 8 correspondence x layer families
                                                #   (no lambda) -> gene_calls_D1full_noLambda.csv: 892 comparisons, 113 at q < 0.05;
                                                #   copy to manuscript_figures/fig23_assets/ for make_fig5_repositioning_v2.py
Rscript test_local.R Igfbp5                     # one gene through the identical code path, on a laptop
```

`01_run_chunk.R`, `02_merge.R` and `03_meta.R` are the earlier layer-pooled sweep (BH within
correspondence × route, blank λ applied). They are kept because `02_merge.R` writes
`gene_spatial_lineage_slices.csv`, the section table that the simulation engine reads (Stage 4);
their `contrasts.csv` is not the table reported in the paper.

The statistic and its analytic null SE (`SE0 = sd(r) · sqrt(Σn² / (Σn)² − 1/N)`) live in
`delta_core.R`; the eligibility gate (≥ 30 cells and `0 < SE0 < 0.05`, `GATE_FREEZE_20260913.md`)
and the routing rules in `inclusion_criteria.R`. The four atlas-supported lineages are declared in
`delta_core.R`. SLURM wrappers (`submit_all.sh`, `*.sbatch`) chain the steps on a cluster.

For the skin atlas (one section per condition) `Skin_MC903_MERFISH/code/skin_application_report.Rmd`
computes the within-section δ and its standardised score without cross-section pooling.

---


### Stage 4 — simulation evaluation

`gene_spatial_cluster/simulation/pilot/engine.R` is the generator and the implementation of the
test and of the baselines (`pilot/LOCKED_METHODS.md` is their frozen specification). `v16/` holds the
experiments reported in the paper; `v17_revision/` the bias/coverage study and the absolute-centroid
baseline. The operative design (10 sections per condition, drawn as four section groups of 2/3/3/2 sections
that the engine labels "mice" — the Figure 4 experiments use no between-group effect; M = 64, NB size
1, τθ ∈ {0, 0.05, 0.15, 0.30}, 200 replicates, 50 FDR panels) is §0 of
`revised_gene_repositioning_simulation_plan_v16.md`.

| Experiment | Script | Output | Figure |
|---|---|---|---|
| τθ → realized I² map; mouse-level variance reference | `v16/02_tau_i2_map_k10.R`, `v16/03_sigma_mouse.R` | `tau_i2_k10.csv`, `sigma_mouse*.csv` | text |
| Type I error, 12 null designs | `v16/run_exp1.R`, `summarise_exp1.R` | `exp1_typeI.csv` | Fig. 4A, 4C |
| Power (gained, reversed) and specificity | `v16/run_batch2.R`, `summarise_batch2.R` | `batch2_power.csv`, `batch2_spec.csv` | Fig. 4B, 4C |
| FDR with 500-gene panels | `v16/run_batch2.R` | `batch2_fdr.csv` | Fig. 4D, supplement |
| REML + Knapp–Hartung vs fixed-effect pooling | `v16/06_power_preview.R` | `ablB_k10.csv` | supplement |
| Robustness (real geometry, abundance, 20 sections, mouse correlation, Fibro 6) | `v16/run_batch3.R`, `summarise_batch3.R` | `batch3_table.csv` | text |
| Absolute-centroid baseline under whole-population movement | `v17_revision/r8_absolute_centroid_baseline.R` | `r8_absolute_centroid_baseline.csv` | Fig. 4C |
| Bias and coverage of Δ̂ | `v17_revision/r3_bias_coverage.R` | `r3_bias_coverage.csv` | text |

The substrate is the real Day 0 Fibro 2 geometry, read by `engine.R::load_substrate()` from the
stage-3 output `gene_spatial_cluster/gene_spatial_lineage_slices.csv` (section → mouse map) and the
MU coordinate files of stage 1.

```bash
cd gene_spatial_cluster/simulation/v16
R_N=5 G0=5 Rscript run_exp1.R 1            # smoke run: cell 1 of 12, 5 replicates x 5 null genes
sbatch exp1.sbatch; sbatch batch2.sbatch; sbatch batch3.sbatch    # full runs (`v16/README_cluster.md`)
Rscript -e 'rmarkdown::render("publication_simulation_report_v16.Rmd")'   # tables from the result CSVs
```

---


### Figures

Run the Level 2 commands; then `python3 manuscript_figures/make_figure_assets.py` is the only extra step
(it rebuilds the figure-input tables from the stage outputs and the atlas before the figure scripts run,
and must follow `make_fig2_correspondence.py`, whose pooled matrix it subsets).
`manuscript_figures/abundance_sensitivity_v2.py` computes the abundance correlations quoted in Section 2.3.

---

## Notes

* Scripts resolve their inputs relative to this folder; `COORD_DIR`, `ALL_SLICES`, `ADATA`, `BENCH` and
  `SIM_DIR` environment variables override the defaults where a script documents them.
* The cluster scripts (`*.sbatch`, `config.sh`, `setup_envs.sh`) target a SLURM system with conda;
  partition names and proxy settings are at the top of each file and must be edited for another cluster.
* Runtimes: a Day 0 × Day 9 correspondence pass over 210 pairs takes minutes on a laptop; the
  gene-level chunk array about an hour per chunk on one core; a simulation cell of 200 replicates
  10–30 min for the proposed test alone and 2–4 h when the related-hypothesis methods implemented in
  `engine.R` (GAM-fs, DESpace2) are scored as well.
* `manuscript_figures/archive/` holds an earlier draft figure script that is not used.
