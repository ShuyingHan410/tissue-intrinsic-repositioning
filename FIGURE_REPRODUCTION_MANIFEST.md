# Figure reproduction manifest

Every figure of the submitted manuscript and supplement, with the script that produces it, its inputs,
where those inputs come from, the output file, the compute needed, and the verification performed on
2026-11-03. "Checkpoint" = the release package of curated coordinates and result tables (see
`reproducibility/README.md`, Level 2). Verification = the script was run in a fresh clone of this
repository with only the checkpoint files added, and its PNG output compared byte-for-byte with the
submitted file. PDF outputs carry an embedded creation timestamp (matplotlib / cairo) and therefore
differ from the submitted PDFs in metadata only.

| Figure | Script (`reproducibility/manuscript_figures/`) | Inputs | Source of inputs | Output | Compute | Status |
|---|---|---|---|---|---|---|
| Fig. 1 framework | `make_fig1_v11.py "$PWD"` | coordinates of two example sections (`output_threelayer/per_slice_morphogam/`), `fig1_assets/igfbp5_cells.csv`, `fig1_assets/pooled_sim_sub.csv` | checkpoint (coordinates, figure_inputs); regenerable by `make_figure_assets.py` from the atlas | `fig1_framework_v11.pdf/.png` | laptop, < 1 min | PNG byte-identical |
| Fig. 2 correspondence | `make_fig2_correspondence.py` | all 93 coordinate files; the 210 per-pair matrices `benchmark_out/day0_day9/morphogam/partials/` | checkpoint (coordinates, benchmark); Level 3: `script/day0_day9/morphogam_recovery.Rmd` | `fig2_correspondence.pdf/.png`; also `fig23_assets/fig2A_pooled_similarity.csv` | laptop, < 1 min | PNG byte-identical; pooled matrix CSV identical |
| Fig. 3 benchmark | `make_fig3_benchmark_v2.py` | `fig23_assets/per_pair_norm_rank_v2.csv`, `benchmark_out/day9_type_counts_full.csv`, `benchmark_out/day9_fibro_totals.csv` | checkpoint (figure_inputs); regenerable by `make_figure_assets.py` from `per_pair_scores_all_methods.csv` (Level 3: cluster runs of moscot, PASTE2, SLAT, STalign + `rescore_all_methods.R`) and the section tables | `fig3_benchmark_v2.pdf/.png` | laptop, < 1 min | PNG byte-identical |
| Fig. 4 simulation | `make_fig4_simulation.R` | `simulation/v16/{exp1_typeI,batch2_power,batch2_spec,batch2_fdr}.csv`, `v17_revision/r8_absolute_centroid_baseline.csv` | checkpoint (simulation); Level 3: `v16/run_exp1.R`, `run_batch2.R`, `v17_revision/r8_*.R` (cluster, days) | `fig4_simulation.pdf/.png`, `fig4C_specificity_values.csv` | laptop, < 1 min | PNG byte-identical |
| Fig. 5 gene repositioning | `make_fig5_repositioning_v2.py` | `fig23_assets/gene_calls_D1full_noLambda.csv`, `fig23_assets/gene_cells_fig5.csv`, `fig1_assets/igfbp5_cells.csv`, coordinates of two sections | checkpoint; `gene_calls` regenerable by `gene_spatial_cluster/04_figure5_calls.R` from `contrastsC.csv` (verified identical), cells by `make_figure_assets.py` | `fig5_repositioning_v2.pdf/.png` | laptop, < 1 min | PNG byte-identical |
| Supp. FDR and ablation | `make_fig3_gene_simulation.R` | `simulation/v16/{batch2_fdr,ablB_k10,exp1_typeI,batch2_power,batch2_spec}.csv` | checkpoint (simulation) | `figS_fdr_ablation.pdf/.png` (+ two `legacy_*` panels not in the submission) | laptop, < 1 min | PNG byte-identical |
| Supp. benchmark sensitivity | `make_figS_benchmark_sensitivity.py` | `benchmark_out/day0_day9/per_pair_scores_all_methods.csv`, `atlas_sections.csv`, coordinates (section list) | checkpoint | `cell_benchmark_sensitivity.pdf/.png`, two CSVs | laptop, < 1 min | PNG byte-identical; CSVs identical |
| Supp. skin application | `make_figS_skin_validation.R` | values transcribed from `Skin_MC903_MERFISH/code/skin_main_results.Rmd` (in the script) | none | `skin_validation_supplement.pdf/.png`, two CSVs | laptop, seconds | PNG byte-identical |
| Supp. animal-level sensitivity (two figures in the current supplement draft) | not in this repository | — | `recomb2027_revision/` (authors' working folder) | — | — | pending the authors' decision on whether the analysis stays in the manuscript |

## Figure-input tables (derived; shipped in the checkpoint and regenerable)

| File | Generator | Validation against the shipped file |
|---|---|---|
| `all_slices/<Slice_ID>.csv` (52), `atlas_sections.csv` | `reproducibility/export_colon_sections.py adata.h5ad` | all 52 tables identical (columns, row order, identifiers, values) |
| `fig1_assets/igfbp5_cells.csv` | `make_figure_assets.py` (atlas `raw.X`) | byte-identical |
| `fig1_assets/pooled_sim_sub.csv` | `make_figure_assets.py` (subset of `fig2A_pooled_similarity.csv`) | byte-identical |
| `fig23_assets/per_pair_norm_rank_v2.csv` | `make_figure_assets.py` (from `per_pair_scores_all_methods.csv`, raw scoring, 590 triples) | byte-identical |
| `fig23_assets/gene_cells_fig5.csv` | `make_figure_assets.py` (coordinates + atlas `raw.X`) | identical content (row order differs) |
| `fig23_assets/gene_calls_D1full_noLambda.csv` | `gene_spatial_cluster/04_figure5_calls.R contrastsC.csv` | identical to machine precision (max |diff| 4e-16), classes identical |
| `fig23_assets/fig2A_pooled_similarity.csv` | `make_fig2_correspondence.py` | byte-identical |
| `benchmark_out/day9_type_counts_full.csv`, `day9_fibro_totals.csv` | `make_figure_assets.py` (section tables) | identical content (quoting differs) |

## Inputs that cannot be regenerated automatically

The 93 coordinate files: 74 are hand-traced MorphoGAM fits (`CurveFinderInteractive`), 13 automatic
and 6 automatic with `knn = 50` (`coordinates_MANIFEST.csv` in the release). `reproducibility/three_layer.Rmd`
reproduces the automatic rungs; the hand-traced curves depend on the author's tracing and are distributed
instead. Everything downstream of the coordinates is deterministic given them.
