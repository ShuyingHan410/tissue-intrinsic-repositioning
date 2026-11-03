Skin (MC903 / oxazolone dermatitis MERFISH, Imanishi et al. 2025) application.

* `data/` (not tracked): the Mendeley Seurat object and the per-cell / count tables written by `code/extract_cells.R` and `code/export_counts.R`.
* `all_slices/` (not tracked): per-strip inputs in the colon schema, written by `code/export_slices.R` (layers from `code/skin_layers.R`).
* `output/morphogam_auto/`: plain MorphoGAM fits per strip x compartment (`code/run_plain_morphogam.R`, then `code/harmonise_orientation.R`). The coordinate files are shipped as `data/skin_tissue_intrinsic_coordinates.tar.gz`; `_AUDIT_plain.csv` and `_ORIENTATION_log.csv` are tracked.
* `code/skin_application_report.Rmd`: the consolidated analysis (correspondence of normal fibroblasts to inflammation-associated fibroblasts, and per-section gene-level delta); it sources the colon engine `script/common/threestep_engine.R` unchanged. `code/overview.Rmd` documents the depth-compartment definitions.
