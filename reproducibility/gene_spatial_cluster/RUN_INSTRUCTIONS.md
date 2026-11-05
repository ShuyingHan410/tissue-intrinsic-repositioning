# How to run the gene spatial-lineage analysis

Measures, per gene per lineage, whether disease changed **where** the lineage deploys that
gene along the crypt→serosa axis.

    delta = sum(n*r)/sum(n) - mean(r)      per gene x cell type x slice, within one layer

`n` = molecule count, `r` = MorphoGAM coordinate (harmonised: lumen = 0 → outer = 1).
delta > 0 = toward crypt base/outer, delta < 0 = toward lumen.
Units = fraction of layer thickness (MU ≈ 266 µm, so 0.025 ≈ 6.7 µm ≈ one cell diameter).

---

## 0. Where things live

| | path |
|---|---|
| this folder (cluster) | `/n/home02/shuyhan/Rong_lab/mophogam/gene_spatial_cluster` |
| expression | `/n/home02/shuyhan/Rong_lab/samples/adata.h5ad` |
| coordinates | `/n/home02/shuyhan/Rong_lab/mophogam/output_threelayer/per_slice_morphogam` |

Both are already set in `config.sh`.

---

## 1. Upload the code (from the Mac)

    cd "/Users/hanshuying/Desktop/PHchan/Lab/Ma lab/benchmark"
    rsync -av --exclude='__pycache__' --exclude='results/local_test_*' \
      gene_spatial_cluster \
      shuyhan@login.rc.fas.harvard.edu:/n/home02/shuyhan/Rong_lab/mophogam/

Coordinates only need re-uploading if the cluster copy predates the 2026-08-10 harmonisation
(`99_check_paths.sh` will tell you):

    rsync -av --delete output_threelayer/per_slice_morphogam \
      shuyhan@login.rc.fas.harvard.edu:/n/home02/shuyhan/Rong_lab/mophogam/output_threelayer/

---

## 2. Set up the environment (once per login session)

    cd /n/home02/shuyhan/Rong_lab/mophogam/gene_spatial_cluster
    module load R                       # R is a module; python is conda (base)

`config.sh` sets `ENV_SETUP`, which every batch job runs. It must match how *you* get python
and R interactively. Verify it works the way a compute node sees it (no interactive shell):

    env -i bash -c 'source config.sh; eval "$ENV_SETUP"; which python3; which Rscript'

Both must print a path. If `Rscript` is missing, pin the exact module version:

    module load R && module list        # note the exact string, e.g. R/4.5.3-fasrc01
    # then edit ENV_SETUP in config.sh to use that exact name

---

## 3. Pre-flight check (always run this first)

    bash 99_check_paths.sh

Expected — **all OK**:

    OK   python3 / Rscript
    OK   ADATA      ... (9.2G)
    OK   COORD_DIR  ... (31 MU files; expect 31)
    OK   raw/X shape (1348582, 990)   |   OK  48 Blank barcodes
    OK   93 files | 0 still INVERTED (want 0) | 5 ambiguous (expected 5)
    OK   python: h5py numpy pandas    |   OK  R: dplyr readr tidyr
    == ALL CHECKS PASSED

The orientation line is the one that matters most: if any file is still inverted, every delta
comes out sign-flipped. Do not proceed until it reads 0.

Also confirm the partition is real for your account:

    sinfo -s | head
    # if 'shared' is wrong:
    sed -i 's/#SBATCH -p shared/#SBATCH -p YOUR_PARTITION/' 0*.sbatch

---

## 4A. Full run (first time, or if the h5ad changed)

    bash submit_all.sh

Chains three jobs with `--dependency=afterok`:

| job | what | resources | time |
|---|---|---|---|
| `00_extract` | h5ad → `meta.csv.gz` + 33 chunks | 1 core, 16 GB | 1–2 h |
| `01_sweep` | **array 1–33**, 30 genes each, parallel | 1 core, 8 GB | 15–45 min |
| `02_merge` | blank calibration + classification | 1 core, 8 GB | ~10 min |

## 4B. Re-score only (code changed, data did not) ← use this after any fix

The extraction is expensive and unchanged, so keep `meta.csv.gz` and `chunks/`:

    rm -f results/part_*.csv results/pairs_*.csv results/slices_*.csv
    rm -f gene_spatial_lineage_*.csv hits.csv

    j=$(sbatch --parsable 01_sweep.sbatch)
    sbatch --dependency=afterok:$j 02_merge.sbatch
    squeue -u shuyhan

**~45 minutes.** Do NOT use `submit_all.sh` here — it would redo the extraction.

---

## 5. Monitor

    squeue -u shuyhan
    tail -f logs/extract_*.out          # during extraction
    tail -f logs/sweep_1.out            # any one array task

**After extraction finishes, check both:**

    ls chunks/*.csv.gz | wc -l          # must be 33
    zcat meta.csv.gz | head -1          # must contain blank_counts

Without `blank_counts` the expression test cannot run. The extractor asserts on it, but confirm.

---

## 6. Check the output is sane BEFORE interpreting

    grep -A6 "blank null for D1" logs/merge_*.out

**Each of the four lineages must show `n_blank` in the dozens.** If a lineage has 0–2 blanks,
its threshold is meaningless and its hits are not trustworthy.

    grep -A8 "routes:" logs/merge_*.out

`excluded` must be present and the four route columns should sum to ~990 genes per lineage.
If genes are missing entirely, rows are being dropped somewhere.

Validation anchor — the paper's one published positional finding (p.2017) must reproduce:

    Rscript -e 'h <- read.csv("hits.csv")
      subset(h, gene %in% c("Mmp13","Mmp10","Plau") & lineage=="Fibro2->IAF4",
             select=c(gene,design,class,delta_d0,delta_d9,median_shift,consistency))'

Those three should be among the strongest luminal (negative delta) genes in IAF 4. If they are
not, something differs from the locally validated path and the rest of the list is suspect.

---

## 7. Outputs

| file | one row per | use |
|---|---|---|
| `gene_spatial_lineage_results.csv` | gene × lineage × design | everything, incl. `route` / `exclude_reason` for untested genes |
| `gene_spatial_lineage_pairs.csv` | **slice pair** | raw: both deltas + shift + slice/mouse provenance — build any summary you like |
| `gene_spatial_lineage_slices.csv` | gene × lineage × slice | rawest: per-slice delta, cells, molecules |
| `hits.csv` | as results, + calibration | `clears_null`, `designs_agree`, `class`, `rank_ok` |

Key columns: `delta_d0` / `delta_d9` (positions), `median_shift`, **`consistency`** (% of slice
pairs agreeing on the sign; blanks ≈ 56% = the floor), `route`, `marker` (annotation only —
does not filter), `class` (preserved / amplified / gained / lost / daughter_positioned / …).

Pull one lineage:

    res <- read.csv("gene_spatial_lineage_results.csv")
    d1  <- subset(res, lineage == "Fibro2->IAF4" & design == "D1_lineage")

---

## 8. Download to the Mac

    rsync -av \
      shuyhan@login.rc.fas.harvard.edu:/n/home02/shuyhan/Rong_lab/mophogam/gene_spatial_cluster/{gene_spatial_lineage_*.csv,hits.csv,genes_index.csv,logs} \
      "/Users/hanshuying/Desktop/PHchan/Lab/Ma lab/benchmark/gene_spatial_cluster/"

---

## Single-gene test (any time, seconds)

    Rscript test_local.R Wnt2b

Runs one gene through the identical code path and writes
`results/local_test_Wnt2b{,_pairs,_slices}.csv`. Requires
`../benchmark_out/gene_demo/gene_counts_gt.csv`.

---

## Limits to state in any write-up

3 Day 0 vs 5 Day 9 mice → **ranked candidates calibrated on blanks, not an FDR-significant
set**; the ~210 slice pairs come from 31 slices and are not independent. The radial axis is
blind to focal lesions. Targeted 940-gene panel → within-panel ranking only. Cell types were
defined by these same genes (position was never a clustering input, so positional claims are
not circular, but `marker` is reported for transparency).
