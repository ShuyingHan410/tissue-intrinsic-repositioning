#!/usr/bin/env python3
"""Step 0b: per-slice x cell-type detection thresholds from the 48 Blank barcodes.

WHY THIS EXISTS. inclusion_criteria.R's expressed_in() implemented the paper's
2.813-sigma rule as

    thr = mean(b_mean) + 2.813 * sd(b_mean)      # b_mean = one value PER SLICE

which takes the SD *across slices* of the per-slice blank mean. The paper's criterion
(Vysochanskij-Petunin) is the mean and SD *across the 48 blank barcodes* within a
cluster - a different distribution. Because slices vary ~100-fold in size, the
across-slice SD is huge (e.g. IAF 4 @ D9: mean 28.8, sd 32.2 -> thr 119.4), and genuinely
expressed genes were routed to `excluded`: Il1b (13.7x above a single blank barcode),
Cxcl5, and Il11 in the SM lineage.

This script computes, for each (Slice_ID, Tier2, day), the 48 per-barcode totals over
the cells of that type in that slice, and writes mean/sd/threshold. Chunk jobs cannot do
this themselves - each holds only its own 30 gene columns, while the blanks are spread
over several chunks.

    python3 00b_blank_thresholds.py      ->  blank_thresholds.csv
"""
import csv, gzip, os, sys
from collections import defaultdict

OUT   = os.environ.get("OUT", ".")
VP_K  = float(os.environ.get("VP_K", "2.813"))

gi = list(csv.DictReader(open(os.path.join(OUT, "genes_index.csv"))))
blank_by_chunk = defaultdict(list)
for r in gi:
    if str(r["is_blank"]).lower() in ("true", "1"):
        blank_by_chunk[int(r["chunk"])].append(r["gene"])
n_blanks = sum(len(v) for v in blank_by_chunk.values())
if not n_blanks:
    sys.exit("no blank barcodes found in genes_index.csv")
print(f"{n_blanks} blank barcodes across chunks {sorted(blank_by_chunk)}", flush=True)

# cell -> (Tier2, day); keep every cell, thresholds are needed per cell type
cell = {}
with gzip.open(os.path.join(OUT, "meta.csv.gz"), "rt") as fh:
    for r in csv.DictReader(fh):
        cell[(r["Slice_ID"], r["Cell_ID"])] = (r["Tier2"], r["day"])
print(f"{len(cell):,} cells in meta", flush=True)

# (slice, type, day) -> {barcode: total}
tot = defaultdict(lambda: defaultdict(float))
ncell = defaultdict(int)
seen_cell = set()
for ch, names in sorted(blank_by_chunk.items()):
    path = os.path.join(OUT, "chunks", f"chunk_{ch:02d}.csv.gz")
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh):
            k = (r["Slice_ID"], r["Cell_ID"])
            m = cell.get(k)
            if m is None:
                continue
            key = (r["Slice_ID"], m[0], m[1])
            if k not in seen_cell:            # count each cell once, not once per chunk
                ncell[key] += 1
                seen_cell.add(k)
            d = tot[key]
            for g in names:
                d[g] += float(r[g])
    print(f"  chunk {ch:02d}: {len(names)} blanks read", flush=True)

def summarise(totals, ncells):
    """mean/sd across barcodes of each barcode's AVERAGE COUNT PER CELL, then the VP bar.
    Per-cell averaging is what the paper specifies and makes the bar comparable across
    groups of different size; multiply by `cells` to get the equivalent raw-total cutoff."""
    vals = [totals[g] / ncells for g in sorted(totals)]
    n = len(vals)
    if n < 2 or ncells <= 0:
        return None
    mean = sum(vals) / n
    sd = (sum((v - mean) ** 2 for v in vals) / (n - 1)) ** 0.5
    return dict(n_blanks=n, blank_mean_per_cell=round(mean, 8),
                blank_sd_per_cell=round(sd, 8),
                thr_per_cell=round(mean + VP_K * sd, 8),
                thr=round((mean + VP_K * sd) * ncells, 4))   # raw-count cutoff for this group

rows = []
# ---- scope = "cluster": the PAPER's definition -----------------------------------------
# STAR Methods p.e13: "we computed the average expression level for each 'blank' barcode
# across all cells within that cluster, then computed the average and standard deviation of
# these values." ONE bar per cluster, pooled over every cell in it. We pool per
# (cluster x day), since the question here is whether a gene is expressed in a cell type at
# a given timepoint - the paper pooled per (cluster x neighborhood) for its ligand-receptor
# analysis, so this is the same construction on a different grouping.
ctot = defaultdict(lambda: defaultdict(float)); cn = defaultdict(int)
for (sl, t, d), bars in tot.items():
    for g, v in bars.items():
        ctot[(t, d)][g] += v
    cn[(t, d)] += ncell[(sl, t, d)]
for key, bars in ctot.items():
    s = summarise(bars, cn[key])
    if s:
        rows.append(dict(scope="cluster", Slice_ID="", Tier2=key[0], day=key[1],
                         cells=cn[key], **s))
# ---- scope = "slice": per-slice variant, kept for robustness comparison ----------------
for key, bars in tot.items():
    s = summarise(bars, ncell[key])
    if s:
        rows.append(dict(scope="slice", Slice_ID=key[0], Tier2=key[1], day=key[2],
                         cells=ncell[key], **s))

rows.sort(key=lambda r: (r["scope"], r["Tier2"], r["day"], r["Slice_ID"]))
path = os.path.join(OUT, "blank_thresholds.csv")
with open(path, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
nc = sum(1 for r in rows if r["scope"] == "cluster")
print(f"\nblank_thresholds.csv: {nc} cluster rows + {len(rows)-nc} slice rows", flush=True)
print("\ncluster-level bars (the paper's definition), per 1000 cells:")
for r in sorted((r for r in rows if r["scope"] == "cluster"), key=lambda r: -r["cells"])[:12]:
    print(f"  {r['Tier2']:10s} {r['day']:6s} cells={r['cells']:7d}  bar/1k cells = {1000*r['thr_per_cell']:7.1f}")
