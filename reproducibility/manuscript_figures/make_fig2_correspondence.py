"""Figure 2 - tissue-intrinsic positional similarity recovers atlas-supported fibroblast correspondences.
A: pooled Day 0 x Day 9 similarity matrix (presence-aware mean over the 210 section pairs, MorphoGAM partials).
B: one Day 0 and one Day 9 section, origin and two candidate populations highlighted (mucosa, Tier2 labels).
C: the same three populations collapsed onto r, medians dashed."""
import os, csv, glob, re, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.dirname(os.path.abspath(__file__)); A = OUT + "/fig23_assets/"
PART = ROOT + "/benchmark_out/day0_day9/morphogam/partials/"; COORD = ROOT + "/output_threelayer/per_slice_morphogam/"
plt.rcParams.update({"font.family": ["Helvetica", "DejaVu Sans"], "font.size": 7, "pdf.fonttype": 42})
NAVY = "#1f3b73"; GREY = "#555555"; DARK = "#333333"; PURPLE = "#5b3fbf"; AMBER = "#e08a1e"; RED = "#c0392b"
D0_ID = "082421_D0_m6_1_slice_1"; D9_ID = "062221_D9_m3_2_slice_2"
ROWS = ["Fibro 6", "Fibro 13", "Fibro 2", "Fibro 7", "Fibro 4", "Fibro 5", "Fibro 12", "Fibro 1", "Fibro 15", "Pericyte 1", "Pericyte 2", "IAF 2"]
COLS = ["Fibro 6", "Fibro 13", "Fibro 2", "Fibro 7", "Fibro 4", "Fibro 5", "Fibro 12", "Fibro 1", "Fibro 15", "Pericyte 1", "Pericyte 2", "IAF 1", "IAF 2", "IAF 3", "IAF 4", "IAF 5", "FRC"]
GT = [("Fibro 2", "IAF 4"), ("Fibro 6", "IAF 2"), ("Fibro 6", "IAF 3"), ("Fibro 7", "Fibro 4")]
MIN_PAIRS = 10

# ---- A: pooled matrix -------------------------------------------------------------
acc = {}; cnt = {}
for f in sorted(glob.glob(PART + "morphogam_Fibroblast_unit*.csv")):
    rr = list(csv.reader(open(f))); hdr = rr[0][1:]
    for r in rr[1:]:
        for c, v in zip(hdr, r[1:]):
            try: x = float(v)
            except ValueError: continue
            if np.isfinite(x): acc[(r[0], c)] = acc.get((r[0], c), 0.0) + x; cnt[(r[0], c)] = cnt.get((r[0], c), 0) + 1
M = np.full((len(ROWS), len(COLS)), np.nan); N = np.zeros_like(M, dtype=int)
for i, r in enumerate(ROWS):
    for j, c in enumerate(COLS):
        if (r, c) in acc: N[i, j] = cnt[(r, c)]; M[i, j] = acc[(r, c)] / cnt[(r, c)] if cnt[(r, c)] >= MIN_PAIRS else np.nan
with open(A + "fig2A_pooled_similarity.csv", "w") as fh:
    w = csv.writer(fh); w.writerow(["Day0"] + COLS)
    for i, r in enumerate(ROWS): w.writerow([r] + [("" if np.isnan(M[i, j]) else f"{M[i, j]:.3f}") for j in range(len(COLS))])
with open(A + "fig2A_n_pairs.csv", "w") as fh:
    w = csv.writer(fh); w.writerow(["Day0"] + COLS)
    for i, r in enumerate(ROWS): w.writerow([r] + list(N[i]))
for o, t in GT:
    i, j = ROWS.index(o), COLS.index(t); row = M[i]; rk = 1 + np.sum(row[np.isfinite(row)] > row[j])
    print(f"{o} -> {t}: pooled {row[j]:.3f}, within-row rank {rk} of {np.isfinite(row).sum()} (n pairs {N[i, j]})")

# ---- B/C: example sections ----------------------------------------------------------
def load(prefix, sid):
    d = {"x": [], "y": [], "r": [], "t2": [], "L": []}
    for L in ("MU", "SM", "ME"):
        for row in csv.DictReader(open(f"{COORD}{prefix}_all_per_slice_{sid}_{L}_xyt_mapped.csv")):
            d["x"].append(float(row["x"])); d["y"].append(float(row["y"])); d["r"].append(float(row["r_for_jsd"])); d["t2"].append(row["Tier2"]); d["L"].append(L)
    return {k: np.array(v) for k, v in d.items()}
d0 = load("Healthy", D0_ID); d9 = load("DSS9", D9_ID)
def kde(vals, grid=np.linspace(0, 1, 300)):
    vals = np.asarray(vals, float); sd = vals.std(); bw = max(0.9 * min(sd, (np.percentile(vals, 75) - np.percentile(vals, 25)) / 1.34) * len(vals) ** (-0.2) * 1.15, 0.02)
    return grid, np.exp(-0.5 * ((grid[:, None] - vals[None, :]) / bw) ** 2).sum(1) / (len(vals) * bw * np.sqrt(2 * np.pi))

fig = plt.figure(figsize=(6.5, 5.9))
def T(x, y, s, size=6.4, color=GREY, ha="center", va="center", **kw): return fig.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, **kw)
def plabel(x, y, s): fig.text(x, y, s, fontsize=11, fontweight="bold", va="top", ha="left")

# A
plabel(0.01, 0.985, "A"); T(0.42, 0.985, "Pooled anatomical similarity: Day 0 fibroblast populations (rows) × Day 9 candidates (columns)", size=7.2, color=NAVY, fontweight="bold")
ax = fig.add_axes((0.115, 0.445, 0.60, 0.44)); cm = plt.get_cmap("viridis").copy(); cm.set_bad("#e6e6e6")
im = ax.imshow(np.ma.masked_invalid(M), cmap=cm, vmin=0, vmax=0.5, aspect="auto")
ax.set_xticks(range(len(COLS))); ax.set_xticklabels(COLS, rotation=45, ha="left", fontsize=6); ax.xaxis.tick_top()
ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(ROWS, fontsize=6); ax.tick_params(length=0, pad=2)
ax.set_xticks(np.arange(-0.5, len(COLS)), minor=True); ax.set_yticks(np.arange(-0.5, len(ROWS)), minor=True); ax.grid(which="minor", color="white", lw=0.8); ax.tick_params(which="minor", length=0)
for sp in ax.spines.values(): sp.set_visible(False)
for i, r in enumerate(ROWS):
    if r in COLS: ax.plot(COLS.index(r), i, "o", ms=2.6, color="white", mec="none")
for o, t in GT:
    ax.add_patch(Rectangle((COLS.index(t) - 0.5, ROWS.index(o) - 0.5), 1, 1, fc="none", ec=RED, lw=1.5, zorder=5))
ax.set_ylabel("Day 0 population", fontsize=6.6, color=DARK, labelpad=3)
cax = fig.add_axes((0.735, 0.60, 0.016, 0.22)); cb = fig.colorbar(im, cax=cax); cb.set_ticks([0, 0.1, 0.2, 0.3, 0.4, 0.5]); cb.ax.tick_params(labelsize=5.6, length=2); cb.outline.set_visible(False)
T(0.743, 0.835, "anatomical\nsimilarity score", size=5.8, color=DARK, ha="left", va="bottom")
LX = 0.80
fig.add_artist(Rectangle((LX, 0.545), 0.016, 0.016, transform=fig.transFigure, fc="none", ec=RED, lw=1.4)); T(LX + 0.024, 0.553, "atlas-supported\ncorrespondence", size=5.8, ha="left", color=DARK)
T(LX, 0.528, "Fibro 2 → IAF 4\nFibro 6 → IAF 2 / IAF 3\nFibro 7 → Fibro 4", size=5.8, ha="left", va="top", color=DARK, linespacing=1.5)
fig.add_artist(Line2D([LX + 0.008], [0.452], transform=fig.transFigure, marker="o", ms=3, color="white", mec="#999999", mew=0.6, ls="none")); T(LX + 0.024, 0.452, "same label in\nboth conditions", size=5.8, ha="left", color=DARK)
fig.add_artist(Rectangle((LX, 0.412), 0.016, 0.016, transform=fig.transFigure, fc="#e6e6e6", ec="none")); T(LX + 0.024, 0.42, f"fewer than {MIN_PAIRS}\ntestable pairs", size=5.8, ha="left", color=DARK)
# B
plabel(0.01, 0.385, "B"); T(0.34, 0.375, "Origin and candidate populations in one Day 0 and one Day 9 section (mucosa)", size=7.2, color=NAVY, fontweight="bold")
def tissue(ax, d, hl, title):
    ax.scatter(d["x"], d["y"], s=0.25, c="#d9d9d9", lw=0, rasterized=True)
    for lab, col in hl: m = d["t2"] == lab; ax.scatter(d["x"][m], d["y"][m], s=1.6, c=col, lw=0, rasterized=True)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=6.4, color=DARK, pad=2, loc="left")
axB0 = fig.add_axes((0.05, 0.05, 0.29, 0.30)); tissue(axB0, d0, [("Fibro 2", DARK)], "Day 0 section")
axB9 = fig.add_axes((0.35, 0.05, 0.29, 0.30)); tissue(axB9, d9, [("IAF 4", PURPLE), ("Fibro 6", AMBER)], "Day 9 section")
hdl = [Line2D([], [], marker="o", ls="none", ms=4, color=DARK, label="Fibro 2 (Day 0 origin)"), Line2D([], [], marker="o", ls="none", ms=4, color=PURPLE, label="IAF 4 (Day 9, atlas-supported match)"), Line2D([], [], marker="o", ls="none", ms=4, color=AMBER, label="Fibro 6 (Day 9, non-match)")]
fig.legend(handles=hdl, loc="lower left", bbox_to_anchor=(0.055, 0.0), ncol=3, fontsize=5.8, frameon=False, handletextpad=0.3, columnspacing=1.2)

# C
plabel(0.655, 0.385, "C"); T(0.835, 0.375, "The same populations on r", size=7.2, color=NAVY, fontweight="bold")
axC = fig.add_axes((0.70, 0.10, 0.27, 0.245))
mu0 = d0["L"] == "MU"; mu9 = d9["L"] == "MU"
sets = [(d0["r"][mu0 & (d0["t2"] == "Fibro 2")], DARK, "Fibro 2"), (d9["r"][mu9 & (d9["t2"] == "IAF 4")], PURPLE, "IAF 4"), (d9["r"][mu9 & (d9["t2"] == "Fibro 6")], AMBER, "Fibro 6")]
ymax = 0
for vals, col, lab in sets:
    g, y = kde(vals); ymax = max(ymax, y.max()); axC.plot(g, y, color=col, lw=1.2); axC.fill_between(g, 0, y, color=col, alpha=0.12, lw=0)
    axC.axvline(np.median(vals), color=col, lw=0.8, ls="--"); print(f"{lab}: n={len(vals)}, median r={np.median(vals):.3f}")
axC.set_xlim(0, 1); axC.set_ylim(0, ymax * 1.12); axC.set_yticks([]); axC.set_xticks([0, 0.5, 1]); axC.set_xticklabels(["0\nlumen", "0.5", "1\nouter"], fontsize=5.8)
for sp in ("top", "right", "left"): axC.spines[sp].set_visible(False)
axC.set_ylabel("cell density", fontsize=6); axC.set_xlabel("r within the mucosa", fontsize=6, labelpad=1)
axC.text(0.02, ymax * 1.08, "dashed: medians", fontsize=5.4, color=GREY, va="top")
fig.savefig(OUT + "/fig2_correspondence.pdf", dpi=300); fig.savefig(OUT + "/fig2_correspondence.png", dpi=300); print("saved")
