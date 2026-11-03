#!/usr/bin/env python3
"""Figure inputs derived from the stage outputs (run after stages 1-3 and make_fig2_correspondence.py).
    python3 make_figure_assets.py            (run from manuscript_figures/; REPRO_ROOT overrides the reproducibility root)
Writes, next to this script:
  fig1_assets/igfbp5_cells.csv          Igfbp5 raw count of every cell of the Figure 1/5 example sections (atlas raw.X)
  fig1_assets/pooled_sim_sub.csv        Fibro 2 / 6 / 7 rows of the pooled similarity matrix (fig23_assets/fig2A_pooled_similarity.csv)
  fig23_assets/per_pair_norm_rank_v2.csv  per-pair rank of the atlas-supported target, raw scoring, levelled to the 590 triples
                                        our method scored (from benchmark_out/day0_day9/per_pair_scores_all_methods.csv)
  fig23_assets/gene_cells_fig5.csv      mucosal Fibro 2 / Fibro 6 (Day 0) and Fibro 2 / Fibro 6 / IAF 3 / IAF 4 (Day 9) cells with
                                        r and the raw counts of Sfrp1, F3, Igfbp5, Cxcl12, Plau (coordinate files + atlas raw.X)
and in <root>/benchmark_out/: day9_type_counts_full.csv and day9_fibro_totals.csv (fibroblast-lineage Tier2 counts of the
21 curated Day 9 sections, from all_slices/)."""
import os, re, csv, collections, h5py, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.environ.get("REPRO_ROOT", os.path.dirname(HERE)); os.chdir(ROOT)
H5 = os.environ.get("ADATA", "adata.h5ad"); COORD = "output_threelayer/per_slice_morphogam"
os.makedirs(f"{HERE}/fig1_assets", exist_ok=True); os.makedirs(f"{HERE}/fig23_assets", exist_ok=True); os.makedirs("benchmark_out", exist_ok=True)
SEC = [r["Slice_ID"] for r in csv.DictReader(open("atlas_sections.csv"))]          # written by export_colon_sections.py
d0 = sorted(s for s in SEC if "_D0_" in s); d9a = sorted(s for s in SEC if "_D9_" in s)
d9m = sorted({f.split("per_slice_")[1].rsplit("_", 3)[0] for f in os.listdir(COORD) if f.startswith("DSS9") and f.endswith("_MU_xyt_mapped.csv")})
assert len(d0) == 10 and len(d9a) == 25 and len(d9m) == 21, (len(d0), len(d9a), len(d9m))

# ---- atlas obs + raw counts of the genes needed
with h5py.File(H5) as h:
    o = h["obs"]
    def col(k):
        g = o[k]; return g["categories"][:].astype(str)[g["codes"][:]] if isinstance(g, h5py.Group) else g[:].astype(str)
    slice_, cid = col("Slice_ID"), col("Cell_ID"); rv = list(h["raw"]["var"]["_index"][:].astype(str))
    GENES = ["Sfrp1", "F3", "Igfbp5", "Cxcl12", "Plau"]; ji = [rv.index(g) for g in GENES]
    sel = np.isin(slice_, d0 + d9m); rows_idx = np.where(sel)[0]
    X = h["raw"]["X"]; counts = np.empty((len(slice_), len(ji)), dtype=np.int32); B = 50000
    for a in range(0, X.shape[0], B):                      # contiguous row blocks: one sequential pass over raw.X
        counts[a:a + B] = X[a:a + B][:, ji].astype(np.int32)
    counts = counts[rows_idx]
cnt = {(s, c): v for s, c, v in zip(slice_[rows_idx], cid[rows_idx], counts)}

# ---- fig1_assets/igfbp5_cells.csv: every cell of the two example sections
EX = ["062921_D0_m3a_1_slice_1", "100221_D9_m5_2_slice_1"]
with open(f"{HERE}/fig1_assets/igfbp5_cells.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["Slice_ID", "Cell_ID", "Igfbp5"])
    for i in rows_idx:
        if slice_[i] in EX: w.writerow([slice_[i], cid[i], cnt[(slice_[i], cid[i])][2]])

# ---- fig23_assets/gene_cells_fig5.csv: mucosal cells of the populations shown in Figure 5, with r
KEEP = {"D0": {"Fibro 2", "Fibro 6"}, "D9": {"Fibro 2", "Fibro 6", "IAF 3", "IAF 4"}}
with open(f"{HERE}/fig23_assets/gene_cells_fig5.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["Slice_ID", "Cell_ID", "day", "Tier2", "r"] + GENES)
    for day, pref, sids in (("D9", "DSS9", d9m), ("D0", "Healthy", d0)):
        for sid in sids:
            for r in csv.DictReader(open(f"{COORD}/{pref}_all_per_slice_{sid}_MU_xyt_mapped.csv")):
                if r["Tier2"] in KEEP[day] and (sid, r["Cell_ID"]) in cnt:
                    w.writerow([sid, r["Cell_ID"], day, r["Tier2"], r["r_for_jsd"]] + list(cnt[(sid, r["Cell_ID"])]))

# ---- fig23_assets/per_pair_norm_rank_v2.csv
LINK = {"Fibro 2": "Fibro 2 → IAF 4", "Fibro 6": "Fibro 6 → IAF 2/IAF 3", "Fibro 7": "Fibro 7 → Fibro 4"}
METHODS = ["morphogam", "moscot", "paste2", "stalign", "slat"]
def pair(r):
    u = int(re.search(r"unit(\d+)", r["unit"]).group(1)); return (d0[u // 21], d9m[u % 21]) if r["method"] == "morphogam" else (d0[u // 25], d9a[u % 25])
gen = {}
for r in csv.DictReader(open("benchmark_out/day0_day9/per_pair_scores_all_methods.csv")):
    if r["scoring"] == "raw" and r["method"] in METHODS: gen[(r["method"],) + pair(r) + (LINK[r["origin"]],)] = (r["rank"], r["n_cand"], r["norm_rank"])
keys = {k[1:] for k in gen if k[0] == "morphogam"}; assert len(keys) == 590
with open(f"{HERE}/fig23_assets/per_pair_norm_rank_v2.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["method", "day0", "day9", "linkage", "rank", "n_cand", "norm_rank"])
    for m in METHODS:
        for k in sorted(keys): w.writerow([m, *k, *gen[(m,) + k]])

# ---- fig1_assets/pooled_sim_sub.csv (needs fig2A_pooled_similarity.csv from make_fig2_correspondence.py)
P = list(csv.DictReader(open(f"{HERE}/fig23_assets/fig2A_pooled_similarity.csv"))); key0 = list(P[0])[0]
SUB = ["IAF 1", "IAF 2", "IAF 3", "IAF 4", "Fibro 4", "Fibro 5"]
with open(f"{HERE}/fig1_assets/pooled_sim_sub.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n"); w.writerow(["row"] + SUB)
    PR = {r[key0]: r for r in P}
    for row in ("Fibro 2", "Fibro 6", "Fibro 7"): w.writerow([row] + [f"{float(PR[row][c]):.3f}" for c in SUB])

# ---- benchmark_out/day9_type_counts_full.csv and day9_fibro_totals.csv
TYPES = ["Fibro 1", "Fibro 12", "Fibro 13", "Fibro 15", "Fibro 2", "Fibro 4", "Fibro 5", "Fibro 6", "Fibro 7", "IAF 1", "IAF 2", "IAF 3", "IAF 4", "IAF 5"]
with open("benchmark_out/day9_type_counts_full.csv", "w", newline="") as fh, open("benchmark_out/day9_fibro_totals.csv", "w", newline="") as ft:
    w = csv.writer(fh); w.writerow(["slice", "type", "n"]); wt = csv.writer(ft, quoting=csv.QUOTE_NONNUMERIC); wt.writerow(["slice", "fibro_total"])
    for sid in d9m:
        c = collections.Counter(r["Tier2"] for r in csv.DictReader(open(f"all_slices/{sid}.csv")))
        for t in TYPES:
            if c[t]: w.writerow([sid, t, c[t]])
        wt.writerow([sid, sum(c[t] for t in TYPES)])
print("written: fig1_assets/igfbp5_cells.csv, fig1_assets/pooled_sim_sub.csv, fig23_assets/per_pair_norm_rank_v2.csv, fig23_assets/gene_cells_fig5.csv, benchmark_out/day9_type_counts_full.csv, benchmark_out/day9_fibro_totals.csv")
