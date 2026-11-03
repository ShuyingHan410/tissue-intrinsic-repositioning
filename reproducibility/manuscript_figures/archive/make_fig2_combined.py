"""Figure 2 (combined) - population correspondence: similarity matrix, example, benchmark, abundance.
A pooled Day 0 x Day 9 similarity matrix; B/C one Day 0 and one Day 9 section and the three populations on r;
D per-pair normalized rank of the atlas-supported match, five methods; E top-three share; F Day 9 composition."""
import os, csv, glob, re, numpy as np, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.dirname(os.path.abspath(__file__)); A = OUT + "/fig23_assets/"
PART = ROOT + "/benchmark_out/day0_day9/morphogam/partials/"; COORD = ROOT + "/output_threelayer/per_slice_morphogam/"
plt.rcParams.update({"font.family": ["Helvetica", "DejaVu Sans"], "font.size": 7, "pdf.fonttype": 42})
NAVY = "#1f3b73"; GREY = "#555555"; DARK = "#333333"; PURPLE = "#5b3fbf"; AMBER = "#e08a1e"; RED = "#c0392b"; TEAL = "#1b7f79"
D0_ID = "082421_D0_m6_1_slice_1"; D9_ID = "062221_D9_m3_2_slice_2"
ROWS = ["Fibro 6", "Fibro 13", "Fibro 2", "Fibro 7", "Fibro 4", "Fibro 5", "Fibro 12", "Fibro 1", "Fibro 15", "Pericyte 1", "Pericyte 2", "IAF 2"]
COLS = ["Fibro 6", "Fibro 13", "Fibro 2", "Fibro 7", "Fibro 4", "Fibro 5", "Fibro 12", "Fibro 1", "Fibro 15", "Pericyte 1", "Pericyte 2", "IAF 1", "IAF 2", "IAF 3", "IAF 4", "IAF 5", "FRC"]
GT = [("Fibro 2", "IAF 4"), ("Fibro 6", "IAF 2"), ("Fibro 6", "IAF 3"), ("Fibro 7", "Fibro 4")]; MIN_PAIRS = 10
METHODS = ["morphogam", "moscot", "paste2", "cast", "stalign"]; MLAB = {"morphogam": "Ours", "moscot": "moscot", "paste2": "PASTE2", "cast": "CAST", "stalign": "STalign"}
LINKS = ["Fibro 6 → IAF 2/IAF 3", "Fibro 7 → Fibro 4", "Fibro 2 → IAF 4"]   # most to least abundant target
TCOL = {"IAF 4": "#9cc3f0", "Fibro 4": "#3a7fd0", "IAF 2/IAF 3": "#0d366b"}; LINK_T = {"Fibro 2 → IAF 4": "IAF 4", "Fibro 6 → IAF 2/IAF 3": "IAF 2/IAF 3", "Fibro 7 → Fibro 4": "Fibro 4"}

fig = plt.figure(figsize=(6.5, 7.3))
def T(x, y, s, size=6.4, color=GREY, ha="center", va="center", **kw): return fig.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, **kw)
def plabel(x, y, s): fig.text(x, y, s, fontsize=11, fontweight="bold", va="top", ha="left")

# ---------------------------------------------------------------- A pooled matrix
acc, cnt = {}, {}
for f in sorted(glob.glob(PART + "morphogam_Fibroblast_unit*.csv")):
    rr = list(csv.reader(open(f))); hdr = rr[0][1:]
    for r in rr[1:]:
        for c, v in zip(hdr, r[1:]):
            try: x = float(v)
            except ValueError: continue
            if np.isfinite(x): acc[(r[0], c)] = acc.get((r[0], c), 0.0) + x; cnt[(r[0], c)] = cnt.get((r[0], c), 0) + 1
M = np.full((len(ROWS), len(COLS)), np.nan)
for i, r in enumerate(ROWS):
    for j, c in enumerate(COLS):
        if cnt.get((r, c), 0) >= MIN_PAIRS: M[i, j] = acc[(r, c)] / cnt[(r, c)]
plabel(0.01, 0.995, "A"); T(0.31, 0.985, "Pooled anatomical similarity, Day 0 populations × Day 9 candidates", size=7.0, color=NAVY, fontweight="bold")
ax = fig.add_axes((0.095, 0.715, 0.44, 0.245)); cm = plt.get_cmap("viridis").copy(); cm.set_bad("#e6e6e6")
im = ax.imshow(np.ma.masked_invalid(M), cmap=cm, vmin=0, vmax=0.5, aspect="auto")
ax.set_xticks(range(len(COLS))); ax.set_xticklabels(COLS, rotation=50, ha="right", fontsize=5.4); ax.set_yticks(range(len(ROWS))); ax.set_yticklabels(ROWS, fontsize=5.4); ax.tick_params(length=0, pad=1.5)
ax.set_xticks(np.arange(-0.5, len(COLS)), minor=True); ax.set_yticks(np.arange(-0.5, len(ROWS)), minor=True); ax.grid(which="minor", color="white", lw=0.7); ax.tick_params(which="minor", length=0)
for sp in ax.spines.values(): sp.set_visible(False)
for i, r in enumerate(ROWS):
    if r in COLS: ax.plot(COLS.index(r), i, "o", ms=2.2, color="white", mec="none")
for o, t in GT: ax.add_patch(Rectangle((COLS.index(t) - 0.5, ROWS.index(o) - 0.5), 1, 1, fc="none", ec=RED, lw=1.3, zorder=5))
ax.set_ylabel("Day 0 population", fontsize=6, color=DARK, labelpad=2); ax.set_xlabel("Day 9 population", fontsize=6, color=DARK, labelpad=1)
cax = fig.add_axes((0.548, 0.80, 0.012, 0.14)); cb = fig.colorbar(im, cax=cax); cb.set_ticks([0, 0.25, 0.5]); cb.ax.tick_params(labelsize=5.2, length=2); cb.outline.set_visible(False)
T(0.548, 0.948, "similarity\nscore", size=5.4, color=DARK, ha="left", va="bottom")
LX = 0.575
fig.add_artist(Rectangle((LX, 0.768), 0.013, 0.013, transform=fig.transFigure, fc="none", ec=RED, lw=1.2)); T(LX + 0.02, 0.7745, "atlas-supported", size=5.3, ha="left", color=DARK)
T(LX, 0.754, "Fibro 2 → IAF 4\nFibro 6 → IAF 2/IAF 3\nFibro 7 → Fibro 4", size=5.3, ha="left", va="top", color=DARK, linespacing=1.45)
fig.add_artist(Line2D([LX + 0.0065], [0.709], transform=fig.transFigure, marker="o", ms=2.6, color="white", mec="#999999", mew=0.6, ls="none")); T(LX + 0.02, 0.709, "same label", size=5.3, ha="left", color=DARK)
fig.add_artist(Rectangle((LX, 0.686), 0.013, 0.013, transform=fig.transFigure, fc="#e6e6e6", ec="none")); T(LX + 0.02, 0.6925, f"< {MIN_PAIRS} pairs", size=5.3, ha="left", color=DARK)

# ---------------------------------------------------------------- B / C example sections
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
plabel(0.665, 0.995, "B"); T(0.835, 0.985, "One Day 0 and one Day 9 section", size=7.0, color=NAVY, fontweight="bold")
def tissue(ax, d, hl, title):
    ax.scatter(d["x"], d["y"], s=0.15, c="#d9d9d9", lw=0, rasterized=True)
    for lab, col in hl: m = d["t2"] == lab; ax.scatter(d["x"][m], d["y"][m], s=1.0, c=col, lw=0, rasterized=True)
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(title, fontsize=5.8, color=DARK, pad=1, loc="center")
tissue(fig.add_axes((0.685, 0.835, 0.15, 0.135)), d0, [("Fibro 2", DARK)], "Day 0")
tissue(fig.add_axes((0.84, 0.835, 0.15, 0.135)), d9, [("IAF 4", PURPLE), ("Fibro 6", AMBER)], "Day 9")
plabel(0.665, 0.825, "C"); T(0.845, 0.815, "The same populations on r (mucosa)", size=7.0, color=NAVY, fontweight="bold")
axC = fig.add_axes((0.715, 0.695, 0.265, 0.11))
mu0 = d0["L"] == "MU"; mu9 = d9["L"] == "MU"
sets = [(d0["r"][mu0 & (d0["t2"] == "Fibro 2")], DARK, "Fibro 2 (Day 0 origin)"), (d9["r"][mu9 & (d9["t2"] == "IAF 4")], PURPLE, "IAF 4 (Day 9, atlas match)"), (d9["r"][mu9 & (d9["t2"] == "Fibro 6")], AMBER, "Fibro 6 (Day 9, non-match)")]
ymax = 0
for vals, col, lab in sets:
    g, y = kde(vals); ymax = max(ymax, y.max()); axC.plot(g, y, color=col, lw=1.1); axC.fill_between(g, 0, y, color=col, alpha=0.12, lw=0)
for vals, col, lab in sets: axC.plot([np.median(vals)] * 2, [0, ymax * 1.12], color=col, lw=0.7, ls="--")
for k, (vals, col, lab) in enumerate(sets): axC.text(0.015, ymax * (1.78 - 0.21 * k), lab, fontsize=5.0, color=col, va="top", fontweight="bold")
axC.set_xlim(0, 1); axC.set_ylim(0, ymax * 1.8); axC.set_yticks([]); axC.set_xticks([0, 0.5, 1]); axC.set_xticklabels(["0 (lumen)", "0.5", "1 (outer)"], fontsize=5.2); axC.tick_params(length=2, pad=1)
for sp in ("top", "right", "left"): axC.spines[sp].set_visible(False)
axC.set_ylabel("cell density", fontsize=5.4, labelpad=1); axC.set_xlabel("r within the mucosa  (dashed: medians)", fontsize=5.4, labelpad=1)

# ---------------------------------------------------------------- D normalized ranks
rows = list(csv.DictReader(open(A + "per_pair_norm_rank.csv")))
nr = collections.defaultdict(list); rk = collections.defaultdict(list)
for r in rows: nr[(r["method"], r["linkage"])].append(float(r["norm_rank"])); rk[(r["method"], r["linkage"])].append(float(r["rank"]))
plabel(0.01, 0.645, "D"); T(0.48, 0.636, "Rank of the atlas-supported Day 9 match, per Day 0 × Day 9 section pair", size=7.0, color=NAVY, fontweight="bold")
rng = np.random.default_rng(1)
for k, L in enumerate(LINKS):
    ax = fig.add_axes((0.075 + k * 0.31, 0.385, 0.27, 0.21))
    data = [nr[(m, L)] for m in METHODS]
    bp = ax.boxplot(data, positions=range(len(METHODS)), widths=0.58, showfliers=False, patch_artist=True, medianprops=dict(color=DARK, lw=1.0), whiskerprops=dict(color="#777777", lw=0.7), capprops=dict(color="#777777", lw=0.7), boxprops=dict(lw=0.7, ec="#777777"))
    for b, m in zip(bp["boxes"], METHODS): b.set_facecolor("#d7ecea" if m == "morphogam" else "#ececec")
    for i, (m, v) in enumerate(zip(METHODS, data)):
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=1.8, color=TEAL if m == "morphogam" else "#8c8c8c", alpha=0.35, lw=0, rasterized=True, zorder=2)
        ax.plot(i, np.mean(v), marker="D", ms=3.2, mfc="white", mec=DARK, mew=0.8, zorder=4)
    ax.set_ylim(1.04, -0.06); ax.set_yticks([0, 0.5, 1]); ax.tick_params(labelsize=5.6, length=2)
    ax.set_xticks(range(len(METHODS))); ax.set_xticklabels([MLAB[m] for m in METHODS], rotation=35, ha="right", fontsize=5.8)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    tc = TCOL[LINK_T[L]]; ax.add_patch(Rectangle((0, 1.02), 1, 0.11, transform=ax.transAxes, fc=tc, ec="none", clip_on=False))
    ax.text(0.5, 1.075, L, transform=ax.transAxes, fontsize=6.4, color="white" if tc != "#9cc3f0" else DARK, ha="center", va="center", fontweight="bold")
    ax.text(0.98, 0.02, f"n = {len(data[0])} pairs" + ("" if len(set(len(d) for d in data)) == 1 else f" (others {len(data[1])})"), transform=ax.transAxes, fontsize=5.0, color=GREY, ha="right", va="bottom")
    if k == 0: ax.set_ylabel("normalized rank\n(0 = ranked first)", fontsize=5.8)
fig.legend(handles=[Line2D([], [], marker="D", ls="none", ms=3.6, mfc="white", mec=DARK, label="mean"), Line2D([], [], color=DARK, lw=1.0, label="median")], loc="upper right", bbox_to_anchor=(0.99, 0.65), fontsize=5.4, frameon=False, ncol=2, handletextpad=0.4, columnspacing=1.0)

# ---------------------------------------------------------------- E top-three share
plabel(0.01, 0.30, "E"); T(0.27, 0.29, "Match within the top three candidates", size=7.0, color=NAVY, fontweight="bold")
ax = fig.add_axes((0.085, 0.105, 0.37, 0.155)); w = 0.15; xs = np.arange(len(LINKS) + 1)
pooled = {m: 100 * np.mean(np.concatenate([np.array(rk[(m, L)]) for L in LINKS]) <= 3) for m in METHODS}
for i, m in enumerate(METHODS):
    vals = [100 * np.mean(np.array(rk[(m, L)]) <= 3) for L in LINKS] + [pooled[m]]
    ax.bar(xs + (i - 2) * w, vals, w, color=TEAL if m == "morphogam" else "#bdbdbd", ec="white", lw=0.3)
    if m == "morphogam":
        for x, v in zip(xs, vals): ax.text(x + (i - 2) * w, v + 2, f"{v:.0f}", fontsize=4.8, ha="center", color=TEAL)
ax.set_xticks(xs); ax.set_xticklabels(["Fibro 6 →\nIAF 2/IAF 3", "Fibro 7 →\nFibro 4", "Fibro 2 →\nIAF 4", "all pairs\npooled"], fontsize=5.4)
ax.set_ylim(0, 140); ax.set_yticks([0, 50, 100]); ax.tick_params(labelsize=5.6, length=2); ax.set_ylabel("% of pairs", fontsize=5.8)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.legend(handles=[Rectangle((0, 0), 1, 1, fc=TEAL, label="Ours"), Rectangle((0, 0), 1, 1, fc="#bdbdbd", label="moscot / PASTE2 / CAST / STalign (in that order)")], fontsize=5.0, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.02), handlelength=1.0, ncol=1)
print("pooled top-3:", {MLAB[m]: round(v, 1) for m, v in pooled.items()})

# ---------------------------------------------------------------- F Day 9 composition
plabel(0.50, 0.30, "F"); T(0.76, 0.29, "Fibroblast-lineage composition of the 21 Day 9 sections", size=7.0, color=NAVY, fontweight="bold")
cnt9 = list(csv.DictReader(open(ROOT + "/benchmark_out/day9_type_counts.csv"))); tot = {r["slice"]: int(float(r["fibro_total"])) for r in csv.DictReader(open(ROOT + "/benchmark_out/day9_fibro_totals.csv"))}
by = collections.defaultdict(lambda: collections.Counter())
for r in cnt9: by[r["slice"]][r["type"]] += int(float(r["n"]))
sl = sorted(tot, key=lambda s: -(by[s]["IAF 2"] + by[s]["IAF 3"]) / tot[s])
comp = {s: {"IAF 2/IAF 3": by[s]["IAF 2"] + by[s]["IAF 3"], "IAF 4": by[s]["IAF 4"], "Fibro 4": by[s]["Fibro 4"]} for s in sl}
for s in sl: comp[s]["other fibroblast lineage"] = tot[s] - sum(comp[s].values())
CC = {"IAF 2/IAF 3": TCOL["IAF 2/IAF 3"], "IAF 4": TCOL["IAF 4"], "Fibro 4": TCOL["Fibro 4"], "other fibroblast lineage": "#dedede"}
ax = fig.add_axes((0.565, 0.105, 0.41, 0.155)); bottom = np.zeros(len(sl))
for key in CC:
    vals = np.array([100 * comp[s][key] / tot[s] for s in sl]); ax.bar(range(len(sl)), vals, 0.78, bottom=bottom, color=CC[key], lw=0); bottom += vals
for i, s in enumerate(sl): ax.text(i, 101, f"{tot[s]/1000:.1f}k", fontsize=3.9, ha="center", va="bottom", color=GREY, rotation=90)
short = lambda s: s.replace("_D9_", " ").replace("_slice_", " s").replace("_", "-")
ax.set_xticks(range(len(sl))); ax.set_xticklabels([short(s) + (" *" if by[s]["IAF 4"] < 30 else "") for s in sl], rotation=90, fontsize=4.1); ax.set_xlim(-0.6, len(sl) - 0.4)
ax.set_ylim(0, 170); ax.set_yticks([0, 50, 100]); ax.tick_params(labelsize=5.6, length=2, pad=1); ax.set_ylabel("% of fibroblast-\nlineage cells", fontsize=5.8)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.legend(handles=[Rectangle((0, 0), 1, 1, fc=CC[k], label=k) for k in CC], fontsize=4.9, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.02), ncol=2, handlelength=1.0, columnspacing=0.8)
fig.savefig(OUT + "/fig2_correspondence_combined.pdf", dpi=300); fig.savefig(OUT + "/fig2_correspondence_combined.png", dpi=300); print("saved")
