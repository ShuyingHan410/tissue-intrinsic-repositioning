# Gene spatial-lineage sweep — FASRC

Measures, per gene x lineage, whether disease changed WHERE the lineage deploys that gene
along the crypt->serosa axis.

    delta = sum(n*r)/sum(n) - mean(r)     per gene x cell type x slice, within one layer

`n` = molecule count, `r` = MorphoGAM coordinate (harmonised 2026-08-10 to lumen=0 -> outer=1).
Sign: delta > 0 = toward crypt base/outer; delta < 0 = toward lumen.
Units: **fraction of layer thickness**, not microns. `r` is normalised 0->1 within each slice's
own layer, so delta is dimensionless and comparable across slices of differing thickness.
(Do not convert to microns or cell diameters with a single global factor: measured MU thickness
ranges 164-962 um across the 31 slices, median 346. Slice thickness does not predict delta —
median |Spearman rho| = 0.16, significant for 1 of 14 top genes — so the normalisation is sound,
but any absolute-distance conversion is not.)

## Analysis reported in the paper (stratified by layer)

    sbatch 01c_stratified.sbatch        # 01c_stratified_chunk.R per chunk -> resultsC/part_NN.csv
    Rscript 02d_stratified_merge.R      # -> contrastsC.csv (892 D1_full comparisons), lambdaC.csv, blankC_by_layer.csv
    Rscript 04_figure5_calls.R          # -> gene_calls_D1full_noLambda.csv (BH within correspondence x layer; 113 at q < 0.05)

Every correspondence x layer is its own analysis; nothing is pooled across layers. The scripts below
(`01_run_chunk.R`, `02_merge.R`, `03_meta.R`) are the earlier layer-pooled sweep, kept for the
section table the simulation engine reads; their `contrasts.csv` / `hits.csv` are not the paper's tables.

## Run on FASRC

Folder location assumed: `/n/home02/shuyhan/Rong_lab/mophogam/gene_spatial_cluster`

    ADATA     = /n/home02/shuyhan/Rong_lab/samples/adata.h5ad
    COORD_DIR = /n/home02/shuyhan/Rong_lab/mophogam/output_threelayer/per_slice_morphogam

1. `bash 99_check_paths.sh`   <- verifies inputs, orientation, and deps BEFORE burning a job
2. Adjust `-p shared` in the three `.sbatch` headers if that is not your partition.
3. `bash submit_all.sh`       (runs the precheck, then chains extract -> sweep -> merge)

After the extract job:  `ls chunks/*.csv.gz | wc -l`  should print **33**.

IMPORTANT: `COORD_DIR` must contain the **harmonised** coordinates (lumen=0 -> outer=1,
harmonised 2026-08-10). `99_check_paths.sh` verifies this and fails if any MU file is still
inverted. Upload the local `output_threelayer/per_slice_morphogam` if the cluster copy is stale.

## Files

| file | role |
|---|---|
| `config.sh` | paths + thresholds (MIN_CELLS, MIN_MOL, MIN_SLICES, MIN_COVERAGE, CHUNK_SIZE) |
| `00_extract_all.py` | h5ad -> `meta.csv.gz` (incl. **blank_counts**) + 33 chunks + `genes_index.csv` |
| `inclusion_criteria.R` | hard exclusions, circularity flag, 2.813-sigma expression test, routing |
| `delta_core.R` | the statistic, slice-pair design, per-pair/per-slice dumps |
| `01_run_chunk.R` | array task: one chunk -> `part_NN.csv`, `pairs_NN.csv`, `slices_NN.csv` |
| `02_merge.R` | bind, blank calibration (per lineage, per route), classify -> `hits.csv` |
| `test_local.R` | `Rscript test_local.R Wnt2b` — one gene through the identical code path |

## Outputs

* `gene_spatial_lineage_results.csv` — one row per gene x lineage x design (nothing dropped;
  `route` / `exclude_reason` say why a gene was not tested)
* `gene_spatial_lineage_pairs.csv` — **one row per slice pair**: both deltas + shift + provenance
* `gene_spatial_lineage_slices.csv` — one row per gene x lineage x slice
* `hits.csv` — calibrated + classified candidates

## Design notes

* **Slice-pair**, not mouse-collapsed: every eligible Day0 slice x every eligible Day9 slice
  (up to 10 x 21 = 210 in MU). Mouse IDs are carried in the pairs table so the dependence can
  be modelled or ignored downstream.
* **Routing**: every gene x lineage is attempted; each gets `D1_full` (shift, both designs),
  `daughter_only` (on at D9 only -> within-daughter position), `parent_only` (lost at D9), or
  `excluded`.
* **Nulls**: the 48 Blank barcodes run through the identical pipeline. D1 rows are judged on
  consistency + |shift|; daughter_only/parent_only on |delta|. Thresholds = blanks' 95th pct,
  per lineage. `02_merge.R` HARD-FAILS if the blank null is empty.
* **`marker`**: annotation only. Flags genes the paper names as defining a lineage's types.
  It does NOT filter anything - every gene goes through the same criteria. The clustering
  used expression and never position, so a defining gene's position inside its own type is
  not circular. (Filtering on it would have dropped Mmp10/Mmp13/Plau.)

## Validated before shipping

| check | source | result |
|---|---|---|
| Wnt2b (base) vs Bmp2 (tip) | Cadinu p.2013 | opposite signs |
| Mmp13/Mmp10/Plau luminal within IAF 4 | Cadinu p.2017 | ranked 1-2-3 of 10, p = 6e-7 |
| epithelial order | textbook | Colonocytes 0.324 < TA 0.467, 26/26 slices |
| blanks | — | ~56% consistency floor |
| Isg15 (level x2.3, no move) | — | 53% -> correctly rejected |

## Known limits

4 Day0 vs 8 Day9 animals -> ranked candidates calibrated on blanks, NOT an FDR-significant set;
210 slice pairs come from 31 slices so they are not independent. The radial axis is blind to
focal lesions. Targeted 940-gene panel -> within-panel ranking only. Cell types were defined
by these same genes (hence the circularity screen). `marker` is a hand-curated list from the
paper's prose, used for annotation only.
