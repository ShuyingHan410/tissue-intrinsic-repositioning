# Reproducibility audit (RECOMB 2027 submission), 2026-11-03

Scope: the public repository `tissue-intrinsic-repositioning` (R package `TICR` + `reproducibility/`)
as of the commit that carries this file. No statistical method, numerical result or manuscript claim
was changed in preparing it; the changes are documentation, path handling, missing generators and
verification. Three levels of reproduction are defined in `reproducibility/README.md`.

## Checks run

| # | Check | How | Result |
|---|---|---|---|
| 1 | Clean installation of the R package | `R CMD INSTALL` into an empty library; `R CMD build` + `R CMD check --no-manual --no-build-vignettes` | **pass** (check status OK) |
| 2 | Unit tests | `testthat::test_dir("tests/testthat")` | **pass** (5 tests, 30 expectations, 0 failures) |
| 3 | Quick-start vignette | `rmarkdown::render("vignettes/quick_start.Rmd")` | **pass** (renders; chunks executed) |
| 4 | README examples | the two quick-start blocks run verbatim in the clean library (`pool_similarity(list(S))`, `test_displacement`) | **pass** (IAF 4 rank 3 of 8 in the single pair; *Igfbp5* reversed, Δ = +0.060) |
| 5 | Figures 1–5 from documented checkpoints | fresh `git clone`; only the checkpoint files added at the paths in Level 2; all five scripts run; PNG compared byte-for-byte with the submitted files | **pass** (5 of 5 byte-identical; PDFs differ in creation timestamp only) |
| 6 | Supplementary figures | same clone: `make_fig3_gene_simulation.R`, `make_figS_benchmark_sensitivity.py`, `make_figS_skin_validation.R` | **pass** (3 of 3 byte-identical). The two animal-level sensitivity figures of the current supplement draft are **not in the repository** (author decision pending) |
| 7 | Figure-input generators | `export_colon_sections.py` vs the 52 section tables used in the paper; `make_figure_assets.py` vs the six shipped figure-input tables; `04_figure5_calls.R` vs the Figure 5 table | **pass** (section tables identical; six of six tables identical in content, four byte-identical, two differing only in row order or quoting; gene calls identical to 4e-16) |
| 8 | Manuscript figure references | `\includegraphics` in `manuscript_revised.tex` / `supplement_revised.tex` mapped to scripts (see `FIGURE_REPRODUCTION_MANIFEST.md`) | **pass** for Figs 1–5 and three supplementary figures; two supplementary figures outside the repository (see 6) |
| 9 | Hard-coded author paths in tracked scripts | `git grep` for home, user and cluster paths | **pass**: none remain; cluster scripts take `BENCH`, `ADATA`, `CONDA_BASE`; stage scripts take `COORD_DIR`, `ALL_SLICES`, `SIM_DIR` |
| 10 | Stage 1 undocumented input | `adata_target_slices_obs_fin.csv` read at the top of `three_layer.Rmd` | **resolved**: the object was never used by the notebook; the line was removed. No replacement file is needed |
| 11 | Level 3 end-to-end run from the atlas | not attempted | **not run** (requires the 6 GB atlas, MorphoGAM, a SLURM cluster for the comparators and the gene sweep, and a step that cannot be automated; see below) |

## What requires author-supplied data or intervention

* **Curated coordinates.** 74 of the 93 coordinate files are hand-traced MorphoGAM fits; 19 are
  automatic. They cannot be rebuilt exactly from the atlas by code. They are prepared as a release
  package (279 MB; manifest with per-file fit provenance and MD5; CC0 source data), **not yet published**;
  Level 2 depends on it. Until it is published, Level 2 is verifiable only by the authors.
* **Comparator benchmark outputs.** moscot, PASTE2, SLAT and STalign were run per section pair on a
  SLURM cluster (6.5 GB of partial outputs, not distributed). The scored per-pair table they feed is in
  the checkpoint release; the harness to recompute them is in `reproducibility/script/`.
* **Gene-level sweep and simulations.** SLURM array jobs (hours per chunk; days for the simulation
  batches). Their result tables are in the checkpoint release.

## Changes made for this audit (no analysis changed)

* Added generators that were missing: `reproducibility/export_colon_sections.py` (atlas → section
  tables), `manuscript_figures/make_figure_assets.py` (six figure-input tables),
  `gene_spatial_cluster/04_figure5_calls.R` (Figure 5 table), `01c_stratified_chunk.R` /
  `02d_stratified_merge.R` (the per-layer gene pipeline behind the paper, previously absent).
* Added the scripts of the submitted figures that were missing or wrong: `make_fig2_correspondence.py`
  (the six-panel draft moved to `manuscript_figures/archive/`), `make_fig4_simulation.R`,
  `make_figS_benchmark_sensitivity.py`, `make_figS_skin_validation.R`; the simulation-figure script's
  non-submitted outputs renamed `legacy_*` so they no longer overwrite Figure 4.
* Removed author-specific paths; `.DS_Store` untracked; `.gitignore` / `.Rbuildignore` restored.
* README: runnable quick start, correct BH family (correspondence × layer), three-level reproduction
  entry point; reproduction guide rewritten with per-step status markers.

## Known limitations stated to the reader

The reproduction guide says explicitly that Level 3 reproduces the automatic coordinate fits but not
the hand-traced ones, that the comparator benchmark and the gene sweep need a cluster, and that the two
animal-level supplementary figures are not yet in the repository.
