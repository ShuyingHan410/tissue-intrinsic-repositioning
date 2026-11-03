"""Figure 5 (v2) - gene repositioning within crypt-associated fibroblast correspondences during DSS colitis.
A significant / tested by layer and correspondence; B signed Delta of the 111 mucosal calls; C the five endpoint
classes; D ONE representative gene (Igfbp5, Fibro 2 -> IAF 4): paired positional profiles pooled over sections;
E Igfbp5 in one Day 0 and one Day 9 section. v2 = single representative gene (was Sfrp1 / Igfbp5 / F3)."""
import os, csv, numpy as np, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
OUT = os.path.dirname(os.path.abspath(__file__)); A = OUT + "/fig23_assets/"; F1 = OUT + "/fig1_assets/"
plt.rcParams.update({"font.family": ["Helvetica", "DejaVu Sans"], "font.size": 7, "pdf.fonttype": 42})
NAVY = "#1f3b73"; GREY = "#555555"; DARK = "#333333"; TEAL = "#1b7f79"; D0C = "#7a7a7a"; D9C = TEAL; BASE = "#1b7f79"; LUM = "#c27c1a"
LAY = ["MU", "SM", "ME"]; LAYN = {"MU": "mucosa", "SM": "submucosa", "ME": "muscularis"}
LIN = ["Fibro2->IAF4", "Fibro6->IAF3", "Fibro6->IAF2", "Fibro7->Fibro4"]; LINN = {l: l.replace("->", " → ").replace("Fibro", "Fibro ").replace("IAF", "IAF ") for l in LIN}
GENE, GLIN, P0, P9 = "Igfbp5", "Fibro2->IAF4", "Fibro 2", "IAF 4"
calls = list(csv.DictReader(open(A + "gene_calls_D1full_noLambda.csv")))
for r in calls:
    for k in ("delta_parent_D0", "delta_daughter_D9", "d1", "q_nolambda"): r[k] = float(r[k])
sig = [r for r in calls if r["q_nolambda"] < 0.05]; print("significant:", len(sig))
W, H = 6.5, 5.15
fig = plt.figure(figsize=(W, H))
def AX(x, y, w, h): return fig.add_axes((x / W, y / H, w / W, h / H))          # positions in inches
def T(x, y, s, size=6.8, color=NAVY, ha="center", va="center", **kw): return fig.text(x / W, y / H, s, fontsize=size, color=color, ha=ha, va=va, fontweight="bold", **kw)
def plabel(x, y, s): fig.text(x / W, y / H, s, fontsize=11, fontweight="bold", va="top", ha="left")
def note(x, y, s, size=4.8, ha="center", va="top", color=GREY): fig.text(x / W, y / H, s, fontsize=size, color=color, ha=ha, va=va)
def clean(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=5.6, length=2)
TOP = H - 0.02

# ---------------------------------------------------------------- A  where: by layer and by correspondence, % of tested
plabel(0.05, TOP, "A"); T(1.30, TOP - 0.07, "Distribution of significant repositioning")
tab = collections.defaultdict(lambda: [0, 0])
for r in calls: tab[(r["lineage"], r["layer"])][1] += 1; tab[(r["lineage"], r["layer"])][0] += r["q_nolambda"] < 0.05
by_layer = {L: [sum(v[0] for k, v in tab.items() if k[1] == L), sum(v[1] for k, v in tab.items() if k[1] == L)] for L in LAY}
by_lin = {l: [sum(v[0] for k, v in tab.items() if k[0] == l), sum(v[1] for k, v in tab.items() if k[0] == l)] for l in LIN}
def hbars(ax, items, names, title, xlab=True):
    ys = np.arange(len(items))[::-1]
    for y, (key, (sg, n)) in zip(ys, items):
        pct = 100 * sg / n; ax.barh(y, pct, 0.62, color="#2a78d6" if pct >= 10 else "#9cc3f0", lw=0)
        ax.text(pct + 1.2, y, f"{sg}/{n}", fontsize=5.2, va="center", color=DARK)
    ax.set_yticks(ys); ax.set_yticklabels([names[k] for k, _ in items], fontsize=5.6); ax.set_xlim(0, 72); ax.set_xticks([0, 25, 50]); ax.tick_params(labelsize=5.4, length=2, pad=1)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    if xlab: ax.set_xlabel("% of tested genes significant", fontsize=5.4, labelpad=1)
    ax.set_title(title, fontsize=5.8, color=DARK, pad=2, loc="left")
hbars(AX(0.70, TOP - 0.85, 1.10, 0.50), [(L, by_layer[L]) for L in LAY], LAYN, "by anatomical layer  (significant / tested)", xlab=False)
hbars(AX(0.70, TOP - 1.78, 1.10, 0.66), [(l, by_lin[l]) for l in LIN], LINN, "by fibroblast correspondence, all layers")

# ---------------------------------------------------------------- B  direction and magnitude (mucosa)
plabel(2.60, TOP, "B"); T(3.68, TOP - 0.07, "Signed Δ of the 111 mucosal calls")
mu = sorted([r for r in sig if r["layer"] == "MU"], key=lambda r: r["d1"])
ax = AX(2.95, TOP - 1.70, 1.45, 1.32); clean(ax); ax.axvline(0, color="#999999", lw=0.7, ls="--")
ys = np.arange(len(mu)); xs = np.array([r["d1"] for r in mu]); cols = [BASE if x > 0 else LUM for x in xs]
ax.scatter(xs, ys, s=5, c=cols, lw=0, zorder=3); ax.set_yticks([]); ax.spines["left"].set_visible(False)
ax.set_xlabel("Δ = δ(Day 9) − δ(Day 0)", fontsize=5.8, labelpad=1); ax.set_xlim(-0.08, 0.1); ax.set_ylim(-3, len(mu) + 2)
npos = int((xs > 0).sum()); nneg = int((xs < 0).sum())
ax.text(0.098, len(mu) * 0.30, f"toward\ncrypt base\nn = {npos}", fontsize=5.4, color=BASE, ha="right", va="center"); ax.text(-0.077, len(mu) * 0.62, f"toward\nlumen\nn = {nneg}", fontsize=5.4, color=LUM, ha="left", va="center")
LBL = {("Sfrp1", "Fibro6->IAF3"): 1, ("F3", "Fibro2->IAF4"): -1, ("Igfbp5", "Fibro2->IAF4"): 1, ("Cxcl12", "Fibro2->IAF4"): 1, ("Plau", "Fibro2->IAF4"): -1}
lab_pts = sorted([(i, r) for i, r in enumerate(mu) if (r["gene"], r["lineage"]) in LBL], key=lambda t: t[0])
slots = {1: [], -1: []}
for i, r in lab_pts:
    sgn = LBL[(r["gene"], r["lineage"])]; y = i
    if slots[sgn] and y - slots[sgn][-1] < 9: y = slots[sgn][-1] + 9
    slots[sgn].append(y); xt = 0.072 if sgn > 0 else -0.052
    is_ex = (r["gene"], r["lineage"]) == (GENE, GLIN)
    if is_ex: ax.scatter([r["d1"]], [i], s=22, facecolor="none", edgecolor=DARK, lw=0.7, zorder=4)
    ax.annotate(r["gene"], (r["d1"], i), (xt, y), fontsize=5.2, style="italic", family="Arial", ha="left" if sgn > 0 else "right", va="center", color=DARK,
                fontweight="bold" if is_ex else "normal", arrowprops=dict(arrowstyle="-", lw=0.4, color="#888888", shrinkA=0, shrinkB=1.5))
ax.set_ylabel("significant tests, ordered by Δ", fontsize=5.6)

# ---------------------------------------------------------------- C  the five endpoint classes, drawn from the data
plabel(4.52, TOP, "C"); T(5.55, TOP - 0.07, "Endpoint classes (δ at Day 0 → Day 9)")
CLS = ["gained", "lost", "reversed", "amplified", "dampened"]
for k, key in enumerate(CLS):
    ax = AX(4.72 + k * 0.355, TOP - 1.52, 0.29, 1.02); members = [r for r in sig if r["pattern"] == key]
    ax.axhline(0, color="#bbbbbb", lw=0.5)
    for r in members:
        col = BASE if r["d1"] > 0 else LUM; ax.plot([0, 1], [r["delta_parent_D0"], r["delta_daughter_D9"]], color=col, lw=0.5, alpha=0.6)
    ax.set_xlim(-0.25, 1.25); ax.set_ylim(-0.11, 0.11); ax.set_xticks([0, 1]); ax.set_xticklabels(["D0", "D9"], fontsize=4.8); ax.set_yticks([-0.1, 0, 0.1] if k == 0 else []); ax.tick_params(labelsize=4.8, length=1.5, pad=1)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    if k > 0: ax.spines["left"].set_visible(False)
    ax.set_title(f"{key}\nn = {len(members)}", fontsize=5.2, color=DARK, pad=2)
    if k == 0: ax.set_ylabel("δ", fontsize=6, style="italic", labelpad=0)
note(5.58, TOP - 1.72, "one line per significant test, coloured by the sign of Δ;\n'positioned' = 95% CI of δ excludes 0; one further test unclassified", size=4.5)

# ---------------------------------------------------------------- D  Igfbp5, pooled over sections
ROW2 = H - 2.08
plabel(0.05, ROW2, "D"); T(0.80, ROW2 - 0.07, GENE, size=6.6, ha="right", style="italic", family="Arial"); T(0.84, ROW2 - 0.07, "in Fibro 2 → IAF 4: transcript position", size=6.4, ha="left")
T(1.45, ROW2 - 0.20, "relative to the host population, pooled over sections", size=6.4)
cells = list(csv.DictReader(open(A + "gene_cells_fig5.csv")))
def kde(vals, w=None, grid=np.linspace(0, 1, 240)):
    vals = np.asarray(vals, float); w = np.ones_like(vals) if w is None else np.asarray(w, float); w = w / w.sum()
    m = np.sum(w * vals); sd = np.sqrt(np.sum(w * (vals - m) ** 2)); neff = 1 / np.sum(w ** 2); bw = max(0.9 * sd * neff ** (-0.2), 0.02)
    return grid, (w[None, :] * np.exp(-0.5 * ((grid[:, None] - vals[None, :]) / bw) ** 2)).sum(1) / (bw * np.sqrt(2 * np.pi)), m
ax = AX(0.50, 1.32, 1.90, 1.38); clean(ax); ax.set_yticks([])
row = [r for r in sig if r["gene"] == GENE and r["lineage"] == GLIN and r["layer"] == "MU"][0]
ymax = 0
for day, pop, col in (("D0", P0, D0C), ("D9", P9, D9C)):
    cc = [c for c in cells if c["day"] == day and c["Tier2"] == pop]; r_ = np.array([float(c["r"]) for c in cc]); w = np.array([float(c[GENE]) for c in cc])
    gx, yh, mh = kde(r_); _, yt, mt = kde(r_[w > 0], w[w > 0]); ymax = max(ymax, yh.max(), yt.max())
    ax.plot(gx, yh, color=col, lw=0.8, ls="--"); ax.plot(gx, yt, color=col, lw=1.3); ax.fill_between(gx, 0, yt, color=col, alpha=0.10, lw=0)
    ax.plot([mh, mh], [0, 0.18], color=col, lw=0.8, ls=":", transform=ax.get_xaxis_transform()); ax.plot([mt, mt], [0, 0.18], color=col, lw=1.0, transform=ax.get_xaxis_transform())
    print(GENE, day, pop, f"cells {len(cc)}, expressing {(w>0).sum()}, pooled delta {mt-mh:+.3f}")
ax.set_xlim(0, 1); ax.set_ylim(0, ymax * 1.40); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "0.5", "1"], fontsize=5.4); ax.set_xlabel("r  (lumen → crypt base)", fontsize=5.6, labelpad=1)
ax.set_ylabel("density along r", fontsize=5.6)
mant, ex_ = f"{row['q_nolambda']:.0e}".split("e"); qs = f"{mant}×10$^{{{int(ex_)}}}$"
ax.text(0.02, 0.97, f"{row['pattern']}, toward {row['direction'].replace('crypt-base', 'crypt base')}\nδ: {row['delta_parent_D0']:+.3f} → {row['delta_daughter_D9']:+.3f}\nΔ = {row['d1']:+.3f},  q = {qs}", transform=ax.transAxes, fontsize=5.2, color=DARK, va="top")
fig.legend(handles=[Line2D([], [], color=D0C, lw=1.3, label="Day 0 (Fibro 2), transcript-weighted"), Line2D([], [], color=D0C, lw=0.8, ls="--", label="Day 0, all cells of the population"),
                    Line2D([], [], color=D9C, lw=1.3, label="Day 9 (IAF 4), transcript-weighted"), Line2D([], [], color=D9C, lw=0.8, ls="--", label="Day 9, all cells of the population")],
           loc="upper left", bbox_to_anchor=(0.30 / W, 1.02 / H), ncol=1, fontsize=5.0, frameon=False, handlelength=1.8, labelspacing=0.35)
note(0.30, 0.34, "ticks at the base: mean positions (solid, transcript-weighted;\ndotted, population); δ is their difference", size=4.8, ha="left")

# ---------------------------------------------------------------- E  Igfbp5 in tissue
plabel(2.60, ROW2, "E"); T(3.32, ROW2 - 0.07, GENE, size=6.6, ha="right", style="italic", family="Arial"); T(3.36, ROW2 - 0.07, "in one Day 0 and one Day 9 section: host population and expressing cells", size=6.4, ha="left")
igf = {(r["Slice_ID"], r["Cell_ID"]): int(float(r["Igfbp5"])) for r in csv.DictReader(open(F1 + "igfbp5_cells.csv"))}
SEC = [("Healthy_062921_D0_m3a_1_slice_1", "062921_D0_m3a_1_slice_1", "Fibro 2", "Day 0 · Fibro 2"), ("DSS9_100221_D9_m5_2_slice_1", "100221_D9_m5_2_slice_1", "IAF 4", "Day 9 · IAF 4")]
vmax = 6
for k, (stem, sid, pop, lab) in enumerate(SEC):
    ax = AX(2.70 + k * 1.76, 0.42, 1.74, 2.15); ax.set_aspect("equal"); ax.axis("off")
    for L in ("ME", "SM", "MU"):
        z = np.load(f"{F1}{stem}_{L}_analysis.npz", allow_pickle=True); ax.scatter(z["x"], z["y"], s=0.12, c="#e3e3e3", lw=0, rasterized=True)
    z = np.load(f"{F1}{stem}_MU_analysis.npz", allow_pickle=True); m = z["Tier3"] == pop; x, y, cid = z["x"][m], z["y"][m], z["Cell_ID"][m]
    cnt = np.array([igf.get((sid, c), 0) for c in cid]); ex = cnt > 0
    ax.scatter(x[~ex], y[~ex], s=1.2, c="#b0b0b0", lw=0, rasterized=True)
    sc = ax.scatter(x[ex], y[ex], s=2.2, c=np.clip(cnt[ex], 0, vmax), cmap="viridis", vmin=1, vmax=vmax, lw=0, rasterized=True, zorder=3)
    ax.set_title(f"{lab}: {m.sum()} cells, {ex.sum()} expressing", fontsize=6, color=DARK, pad=1, y=0.97)
cax = AX(6.24, 1.05, 0.08, 0.90); cb = fig.colorbar(sc, cax=cax); cb.set_ticks([1, 3, 6]); cb.set_ticklabels(["1", "3", "≥6"]); cb.ax.tick_params(labelsize=5.2, length=2); cb.outline.set_visible(False)
fig.text(6.28 / W, 2.00 / H, "Igfbp5\ncounts", fontsize=5.4, color=DARK, ha="center", va="bottom", style="italic", family="Arial")
fig.legend(handles=[Line2D([], [], marker="o", ls="none", ms=3, color="#b0b0b0", label="host population, not expressing"), Line2D([], [], marker="o", ls="none", ms=3, color="#e3e3e3", label="other cells")],
           loc="lower center", bbox_to_anchor=(4.55 / W, 0.02 / H), ncol=2, fontsize=5.0, frameon=False)
fig.savefig(OUT + "/fig5_repositioning_v2.pdf", dpi=300); fig.savefig(OUT + "/fig5_repositioning_v2.png", dpi=300); print("saved")
