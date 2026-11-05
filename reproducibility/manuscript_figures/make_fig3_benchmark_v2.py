"""Figure 3 (v2, 2026-10-31: SLAT added; comparators restricted to exactly the (Day0, Day9, origin) pairs MorphoGAM scored: 200/210/180; raw scoring) - benchmarking population correspondence against alignment-based methods.
A: per-section-pair normalized rank of the atlas-supported Day 9 match, three linkages x five methods.
B: share of testable pairs with the match in the top three candidates.
C: fibroblast-lineage composition of the 21 Day 9 sections, every population shown (recovery targets in blue, others muted)."""
import os, csv, numpy as np, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.dirname(os.path.abspath(__file__)); A = OUT + "/fig23_assets/"
plt.rcParams.update({"font.family": ["Helvetica", "DejaVu Sans"], "font.size": 7, "pdf.fonttype": 42})
NAVY = "#1f3b73"; GREY = "#555555"; DARK = "#333333"; TEAL = "#1b7f79"
METHODS = ["morphogam", "moscot", "slat", "paste2", "stalign"]; MLAB = {"morphogam": "Ours", "moscot": "moscot", "slat": "SLAT", "paste2": "PASTE2", "stalign": "STalign"}
LINKS = ["Fibro 6 → IAF 2/IAF 3", "Fibro 7 → Fibro 4", "Fibro 2 → IAF 4"]  # most to least abundant target: dark → light
# target-population colours shared by the facet headers in A and the composition bars in C (light = rarer target)
TCOL = {"IAF 4": "#9cc3f0", "Fibro 4": "#3a7fd0", "IAF 2/IAF 3": "#0d366b"}; LINK_T = {"Fibro 2 → IAF 4": "IAF 4", "Fibro 6 → IAF 2/IAF 3": "IAF 2/IAF 3", "Fibro 7 → Fibro 4": "Fibro 4"}
MCOL = {m: (TEAL if m == "morphogam" else "#8c8c8c") for m in METHODS}
rows = list(csv.DictReader(open(A + "per_pair_norm_rank_v2.csv")))
nr = collections.defaultdict(list); rk = collections.defaultdict(list)
for r in rows: nr[(r["method"], r["linkage"])].append(float(r["norm_rank"])); rk[(r["method"], r["linkage"])].append(float(r["rank"]))

fig = plt.figure(figsize=(6.5, 5.3))
def T(x, y, s, size=6.4, color=GREY, ha="center", va="center", **kw): return fig.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, **kw)
def plabel(x, y, s): fig.text(x, y, s, fontsize=11, fontweight="bold", va="top", ha="left")

# ---- A: normalized-rank distributions ----------------------------------------------
plabel(0.01, 0.985, "A"); T(0.5, 0.975, "Rank of the atlas-supported Day 9 match, per Day 0 × Day 9 section pair", size=7.4, color=NAVY, fontweight="bold")
rng = np.random.default_rng(1)
for k, L in enumerate(LINKS):
    ax = fig.add_axes((0.075 + k * 0.31, 0.555, 0.27, 0.36))
    data = [nr[(m, L)] for m in METHODS]
    bp = ax.boxplot(data, positions=range(len(METHODS)), widths=0.58, showfliers=False, patch_artist=True, medianprops=dict(color=DARK, lw=1.0), whiskerprops=dict(color="#777777", lw=0.7), capprops=dict(color="#777777", lw=0.7), boxprops=dict(lw=0.7, ec="#777777"))
    for b, m in zip(bp["boxes"], METHODS): b.set_facecolor("#d7ecea" if m == "morphogam" else "#ececec")
    for i, (m, v) in enumerate(zip(METHODS, data)):
        ax.scatter(i + rng.uniform(-0.17, 0.17, len(v)), v, s=2.2, color=MCOL[m], alpha=0.35, lw=0, rasterized=True, zorder=2)
        ax.plot(i, np.mean(v), marker="D", ms=3.6, mfc="white", mec=DARK, mew=0.8, zorder=4)
    ax.set_ylim(1.04, -0.06); ax.set_yticks([0, 0.25, 0.5, 0.75, 1]); ax.tick_params(labelsize=5.8, length=2)
    ax.set_xticks(range(len(METHODS))); ax.set_xticklabels([MLAB[m] for m in METHODS], rotation=35, ha="right", fontsize=6)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    tc = TCOL[LINK_T[L]]; ax.add_patch(Rectangle((0, 1.015), 1, 0.075, transform=ax.transAxes, fc=tc, ec="none", clip_on=False))
    ax.text(0.5, 1.052, L, transform=ax.transAxes, fontsize=6.6, color="white" if tc != "#9cc3f0" else DARK, ha="center", va="center", fontweight="bold")
    ax.text(0.98, 0.02, f"n = {len(data[0])} pairs" + ("" if len(set(len(d) for d in data)) == 1 else f" (ours {len(data[0])}, others {len(data[1])})"), transform=ax.transAxes, fontsize=5.2, color=GREY, ha="right", va="bottom")
    if k == 0: ax.set_ylabel("normalized rank of the match\n(0 = ranked first)", fontsize=6.2)
hdl = [Line2D([], [], marker="D", ls="none", ms=4, mfc="white", mec=DARK, label="mean"), Line2D([], [], color=DARK, lw=1.0, label="median")]
fig.legend(handles=hdl, loc="upper right", bbox_to_anchor=(0.99, 0.995), fontsize=5.6, frameon=False, ncol=2, handletextpad=0.4, columnspacing=1.0)

# ---- B: top-3 share --------------------------------------------------------------------
plabel(0.01, 0.405, "B"); T(0.26, 0.395, "Match within the top three candidates", size=7.2, color=NAVY, fontweight="bold")
ax = fig.add_axes((0.085, 0.145, 0.37, 0.215))
w = 0.15; xs = np.arange(len(LINKS) + 1)
pooled = {m: 100 * np.mean(np.concatenate([np.array(rk[(m, L)]) for L in LINKS]) <= 3) for m in METHODS}
for i, m in enumerate(METHODS):
    vals = [100 * np.mean(np.array(rk[(m, L)]) <= 3) for L in LINKS] + [pooled[m]]
    ax.bar(xs + (i - 2) * w, vals, w, color=TEAL if m == "morphogam" else "#bdbdbd", ec="white", lw=0.3)
    if m == "morphogam":
        for x, v in zip(xs, vals): ax.text(x + (i - 2) * w, v + 1.5, f"{v:.0f}", fontsize=5, ha="center", color=TEAL)
ax.set_xticks(xs); ax.set_xticklabels(["Fibro 6 →\nIAF 2/IAF 3", "Fibro 7 →\nFibro 4", "Fibro 2 →\nIAF 4", "all pairs\npooled"], fontsize=5.8)
ax.set_ylim(0, 140); ax.set_yticks([0, 25, 50, 75, 100]); ax.tick_params(labelsize=5.8, length=2); ax.set_ylabel("% of pairs", fontsize=6.2)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
hdl = [Rectangle((0, 0), 1, 1, fc=TEAL, label="Ours")] + [Rectangle((0, 0), 1, 1, fc="#bdbdbd", label="moscot / SLAT / PASTE2 / STalign (in that order)")]
ax.legend(handles=hdl, fontsize=5.4, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.0), handlelength=1.0, ncol=1)
print("pooled top-3:", {MLAB[m]: round(v, 1) for m, v in pooled.items()})

# ---- C: Day 9 composition (every fibroblast-lineage population; pericytes and FRC are not part of the lineage total) ----
plabel(0.50, 0.405, "C"); T(0.76, 0.395, "Fibroblast-lineage composition of the 21 Day 9 sections", size=7.2, color=NAVY, fontweight="bold")
cnt = list(csv.DictReader(open(ROOT + "/benchmark_out/day9_type_counts_full.csv"))); tot = {r["slice"]: int(float(r["fibro_total"])) for r in csv.DictReader(open(ROOT + "/benchmark_out/day9_fibro_totals.csv"))}
by = collections.defaultdict(lambda: collections.Counter())
for r in cnt: by[r["slice"]][r["type"]] += int(float(r["n"]))
for s in tot: assert sum(by[s].values()) == tot[s], s
sl = sorted(tot, key=lambda s: -(by[s]["IAF 2"] + by[s]["IAF 3"]) / tot[s])
TARGETS = ["IAF 2/IAF 3", "IAF 4", "Fibro 4"]
others = sorted({t for s in sl for t in by[s]} - {"IAF 2", "IAF 3", "IAF 4", "Fibro 4"}, key=lambda t: -sum(by[s][t] / tot[s] for s in sl))
comp = {s: {"IAF 2/IAF 3": by[s]["IAF 2"] + by[s]["IAF 3"], **{t: by[s][t] for t in ["IAF 4", "Fibro 4"] + others}} for s in sl}
# muted, light colours for the non-target populations so the three targets (blues) stay in front
OTHER_COL = ["#d9d9d9", "#e9ddd0", "#cfdbd3", "#e4d7e4", "#c8cfd9", "#e8e3c9", "#d8cabf", "#bccfca", "#f0f0f0", "#b5b5b5", "#dcd3c3", "#cbd7df", "#e2e2e2"]
CC = {"IAF 2/IAF 3": TCOL["IAF 2/IAF 3"], "IAF 4": TCOL["IAF 4"], "Fibro 4": TCOL["Fibro 4"], **{t: OTHER_COL[i % len(OTHER_COL)] for i, t in enumerate(others)}}
ax = fig.add_axes((0.565, 0.145, 0.41, 0.185)); bottom = np.zeros(len(sl))
for key in CC:
    vals = np.array([100 * comp[s][key] / tot[s] for s in sl]); ax.bar(range(len(sl)), vals, 0.78, bottom=bottom, color=CC[key], lw=0.25, ec="white"); bottom += vals
for i, s in enumerate(sl):
    ax.text(i, 101, f"{tot[s]/1000:.1f}k", fontsize=4.2, ha="center", va="bottom", color=GREY, rotation=90)
short = lambda s: s.replace("_D9_", " ").replace("_slice_", " s").replace("_", "-")
ax.set_xticks(range(len(sl))); ax.set_xticklabels([short(s) + (" *" if by[s]["IAF 4"] < 30 else "") for s in sl], rotation=90, fontsize=4.4); ax.set_xlim(-0.6, len(sl) - 0.4)
ax.set_ylim(0, 122); ax.set_yticks([0, 25, 50, 75, 100]); ax.tick_params(labelsize=5.8, length=2, pad=1); ax.set_ylabel("% of fibroblast-lineage cells", fontsize=6.2)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
hdl = [Rectangle((0, 0), 1, 1, fc=CC[k], ec="#bbbbbb", lw=0.3, label=k) for k in CC]
ax.legend(handles=hdl, fontsize=4.6, frameon=False, loc="lower left", bbox_to_anchor=(-0.02, 1.0), ncol=5, handlelength=0.9, handleheight=0.9, columnspacing=0.7, labelspacing=0.25, handletextpad=0.4, borderaxespad=0)
with open(A + "fig3C_day9_composition.csv", "w") as fh:
    w = csv.writer(fh); w.writerow(["slice", "fibro_total"] + list(CC)); [w.writerow([s, tot[s]] + [comp[s][k] for k in CC]) for s in sl]
fig.savefig(OUT + "/fig3_benchmark_v2.pdf", dpi=300); fig.savefig(OUT + "/fig3_benchmark_v2.png", dpi=300); print("saved")
