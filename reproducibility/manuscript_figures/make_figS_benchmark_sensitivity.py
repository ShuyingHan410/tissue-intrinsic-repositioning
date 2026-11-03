#!/usr/bin/env python3
"""Supplementary figure: cell-population benchmark sensitivity (no re-run of any method).
Run from this folder:  python3 make_figS_benchmark_sensitivity.py        (REPRO_ROOT overrides the reproducibility root)
Inputs (relative to the reproducibility root): benchmark_out/day0_day9/per_pair_scores_all_methods.csv (raw / column-normalised /
ties-averaged scoring of the SAME partials, written by script/day0_day9/rescore_all_methods.R), atlas_sections.csv (Mouse_ID per
section, written by export_colon_sections.py) and output_threelayer/per_slice_morphogam/ (the 21 curated Day 9 sections). Comparators are levelled to the exact 590
(origin, Day 0, Day 9) triples our method scored (as in make_fig3_benchmark_v2.py).
Outputs next to this script: cell_benchmark_normalization.csv, cell_benchmark_animalpair.csv, cell_benchmark_sensitivity.pdf/.png"""
import csv, collections, statistics, os, re
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import numpy as np
OUT = os.path.dirname(os.path.abspath(__file__)); os.chdir(os.environ.get("REPRO_ROOT", os.path.dirname(OUT)))
atlas = {r["Slice_ID"]: r for r in csv.DictReader(open("atlas_sections.csv"))}      # written by export_colon_sections.py
d0 = sorted(s for s in atlas if "_D0_" in s); d9a = sorted(s for s in atlas if "_D9_" in s)
d9m = sorted({f.split("per_slice_")[1].rsplit("_", 3)[0] for f in os.listdir("output_threelayer/per_slice_morphogam") if "DSS9" in f})
assert len(d0) == 10 and len(d9a) == 25 and len(d9m) == 21
rows = list(csv.DictReader(open("benchmark_out/day0_day9/per_pair_scores_all_methods.csv")))
METHODS = ["morphogam", "moscot", "slat", "paste2", "stalign"]; LAB = {"morphogam": "Ours", "moscot": "moscot", "slat": "SLAT", "paste2": "PASTE2", "stalign": "STalign"}
SC = ["raw", "column-normalised", "column-normalised, ties averaged"]
def pair(r):
    u = int(re.search(r"unit(\d+)", r["unit"]).group(1))
    return (d0[u // 21], d9m[u % 21]) if r["method"] == "morphogam" else (d0[u // 25], d9a[u % 25])
tab = collections.defaultdict(dict)   # (method, scoring) -> {(d0,d9,origin): (rank, norm)}
for r in rows:
    if r["method"] not in METHODS: continue
    tab[(r["method"], r["scoring"])][pair(r) + (r["origin"],)] = (float(r["rank"]), float(r["norm_rank"]))
keys = set(tab[("morphogam", "raw")]); assert len(keys) == 590
LINK = {"Fibro 2": "Fibro 2 → IAF 4", "Fibro 6": "Fibro 6 → IAF 2/IAF 3", "Fibro 7": "Fibro 7 → Fibro 4"}
def metrics(d):
    rk = [v[0] for v in d.values()]; nr = [v[1] for v in d.values()]
    return dict(n=len(d), top1=100 * np.mean(np.array(rk) <= 1), top3=100 * np.mean(np.array(rk) <= 3), mean_norm_rank=np.mean(nr))
out = []
for m in METHODS:
    for s in SC:
        d = {k: v for k, v in tab[(m, s)].items() if k in keys}; assert len(d) == 590, (m, s, len(d))
        for L in ["all"] + list(LINK):
            dd = d if L == "all" else {k: v for k, v in d.items() if k[2] == L}
            mt = metrics(dd); out.append(dict(method=LAB[m], scoring=s, correspondence="pooled" if L == "all" else LINK[L], **{k: round(v, 3) if k == "mean_norm_rank" else round(v, 1) if k != "n" else v for k, v in mt.items()}))
with open(f"{OUT}/cell_benchmark_normalization.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
# ---- B2: per animal pair ----
ap = []
for m in METHODS:
    d = {k: v for k, v in tab[(m, "raw")].items() if k in keys}
    by = collections.defaultdict(list)
    for (a, b, o), v in d.items(): by[(atlas[a]["Mouse_ID"], atlas[b]["Mouse_ID"])].append(v)
    per = {k: dict(n=len(v), top3=100 * np.mean([x[0] <= 3 for x in v]), top1=100 * np.mean([x[0] <= 1 for x in v]), mnr=np.mean([x[1] for x in v])) for k, v in by.items()}
    for (a, b), v in sorted(per.items()): ap.append(dict(method=LAB[m], mouse_day0=a, mouse_day9=b, n_evaluations=v["n"], top1=round(v["top1"], 1), top3=round(v["top3"], 1), mean_norm_rank=round(v["mnr"], 3)))
    ap.append(dict(method=LAB[m], mouse_day0="MEAN over 15 animal pairs (equal weight)", mouse_day9="", n_evaluations=sum(v["n"] for v in per.values()),
                   top1=round(np.mean([v["top1"] for v in per.values()]), 1), top3=round(np.mean([v["top3"] for v in per.values()]), 1), mean_norm_rank=round(np.mean([v["mnr"] for v in per.values()]), 3)))
    ap.append(dict(method=LAB[m], mouse_day0="RANGE over animal pairs (top3 min-max)", mouse_day9="", n_evaluations="", top1="", top3=f"{min(v['top3'] for v in per.values()):.0f}-{max(v['top3'] for v in per.values()):.0f}", mean_norm_rank=""))
    ap.append(dict(method=LAB[m], mouse_day0="POOLED over 590 evaluations", mouse_day9="", n_evaluations=590, **{k: round(v, 1) if k != "mean_norm_rank" else round(v, 3) for k, v in metrics(d).items() if k != "n"}))
with open(f"{OUT}/cell_benchmark_animalpair.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ap[0].keys())); w.writeheader(); w.writerows(ap)
# ---- figure: A top-3 by scoring; B per-animal-pair top-3 (raw) ----
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.9)); plt.rcParams.update({"font.size": 7})
x = np.arange(len(METHODS)); w = 0.26; cols = ["#1b7f79", "#8fbfbb", "#c9c9c9"]
for i, s in enumerate(SC):
    vals = [next(o["top3"] for o in out if o["method"] == LAB[m] and o["scoring"] == s and o["correspondence"] == "pooled") for m in METHODS]
    ax[0].bar(x + (i - 1) * w, vals, w, color=cols[i], label=s)
ax[0].set_xticks(x); ax[0].set_xticklabels([LAB[m] for m in METHODS], fontsize=7); ax[0].set_ylabel("% of 590 evaluations, match in top 3", fontsize=7); ax[0].set_ylim(0, 100)
ax[0].legend(fontsize=5.5, frameon=False, loc="lower left"); ax[0].set_title("A  Scoring of the comparators' coupling mass", fontsize=7.5, loc="left")
for sp in ("top", "right"): ax[0].spines[sp].set_visible(False)
rng = np.random.default_rng(0)
for i, m in enumerate(METHODS):
    v = [r["top3"] for r in ap if r["method"] == LAB[m] and r["mouse_day0"].startswith("0")]; n = [r["n_evaluations"] for r in ap if r["method"] == LAB[m] and r["mouse_day0"].startswith("0")]
    ax[1].scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=np.array(n) * 0.5, color="#1b7f79" if m == "morphogam" else "#8c8c8c", alpha=0.6, lw=0)
    mean_ap = next(r["top3"] for r in ap if r["method"] == LAB[m] and r["mouse_day0"].startswith("MEAN")); pooled = next(r["top3"] for r in ap if r["method"] == LAB[m] and r["mouse_day0"].startswith("POOLED"))
    ax[1].plot([i - 0.3, i + 0.3], [mean_ap, mean_ap], color="k", lw=1.2); ax[1].plot(i, pooled, marker="D", ms=4, mfc="white", mec="k")
ax[1].set_xticks(x); ax[1].set_xticklabels([LAB[m] for m in METHODS], fontsize=7); ax[1].set_ylim(0, 105); ax[1].set_ylabel("% top-3 within each animal pair", fontsize=7)
ax[1].set_title("B  Per animal pair (3 Day 0 × 5 Day 9 = 15 pairs; raw scoring)", fontsize=7.5, loc="left")
for sp in ("top", "right"): ax[1].spines[sp].set_visible(False)
from matplotlib.lines import Line2D
ax[1].legend(handles=[Line2D([], [], color="k", lw=1.2, label="mean of animal pairs"), Line2D([], [], marker="D", ls="none", mfc="white", mec="k", label="pooled over 590"), Line2D([], [], marker="o", ls="none", color="#8c8c8c", alpha=.6, label="one animal pair (size ∝ evaluations)")], fontsize=5.5, frameon=False, loc="lower left")
fig.tight_layout(); fig.savefig(f"{OUT}/cell_benchmark_sensitivity.pdf"); fig.savefig(f"{OUT}/cell_benchmark_sensitivity.png", dpi=300)
for o in out:
    if o["correspondence"] == "pooled": print(o)
print()
for r in ap:
    if not r["mouse_day0"].startswith("0"): print(r)
