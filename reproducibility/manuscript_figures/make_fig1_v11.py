import sys, os, csv, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Polygon
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe

OUT = sys.argv[1]; A = OUT + "/fig1_assets/"
ROOT = os.path.dirname(OUT.rstrip("/")); COORD = ROOT + "/output_threelayer/per_slice_morphogam/"
plt.rcParams.update({"font.family": ["Helvetica", "DejaVu Sans"], "font.size": 7, "pdf.fonttype": 42})
NAVY="#1f3b73"; GREY="#555555"; RED="#c0392b"; DARK="#333333"; PURPLE="#5b3fbf"; ORANGE="#e08a1e"; TEAL="#1b7f79"
LC = {"MU": "#79b4d6", "SM": "#9b89cf", "ME": "#f0a67a"}; LN = {"MU": "mucosa", "SM": "submucosa", "ME": "muscularis"}
LCM = {"MU": "Blues", "SM": "Purples", "ME": "Oranges"}; LMUTED = {"MU": "#a9cde6", "SM": "#bdb2dc", "ME": "#f4c3a4"}
CURVE = "#4a4a4a"
D0 = ("Healthy", "062921_D0_m3a_1_slice_1"); D9 = ("DSS9", "100221_D9_m5_2_slice_1")
STAGE = ["#8e6cc7", "#e39b3a", "#3f9d6b"]

def load_layer(cond, slide, L):
    cache = f"{A}/{cond}_{slide}_{L}_analysis.npz"
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True); return {k: z[k] for k in z.files}
    num = {k: [] for k in ("x", "y", "r_for_jsd", "f1", "f2", "t")}; cid = []; t3 = []
    for row in csv.DictReader(open(f"{COORD}{cond}_all_per_slice_{slide}_{L}_xyt_mapped.csv")):
        for k in num: num[k].append(float(row[k]))
        cid.append(row["Cell_ID"]); t3.append(row["Tier3"])
    d = {k: np.array(v) for k, v in num.items()}; d["Cell_ID"] = np.array(cid); d["Tier3"] = np.array(t3); np.savez(cache, **d); return d
SEC = {"D0": {L: load_layer(*D0, L) for L in ("MU", "SM", "ME")}, "D9": {L: load_layer(*D9, L) for L in ("MU", "SM", "ME")}}
def orientation_ok(D):
    allxy = np.vstack([np.c_[d["x"], d["y"]] for d in D.values()]); c = allxy.mean(0); out = {}
    for L, d in D.items():
        rad = np.hypot(d["x"]-c[0], d["y"]-c[1]); ang = np.arctan2(d["y"]-c[1], d["x"]-c[0]); bins = np.digitize(ang, np.linspace(-np.pi, np.pi, 37)); cors = []
        for b in np.unique(bins):
            m = bins == b
            if m.sum() > 30 and np.std(rad[m]) > 0 and np.std(d["r_for_jsd"][m]) > 0: cors.append(np.corrcoef(rad[m], d["r_for_jsd"][m])[0, 1])
        out[L] = float(np.median(cors))
    return out
for key in SEC:
    for L, med in orientation_ok(SEC[key]).items():
        if med < -0.3:
            print(f"WARNING: {key} {L} r decreases outward (median corr {med:+.2f}); flipped for DISPLAY ONLY in this figure")
            SEC[key][L]["r_for_jsd"] = 1 - SEC[key][L]["r_for_jsd"]
igf = {(r["Slice_ID"], r["Cell_ID"]): int(float(r["Igfbp5"])) for r in csv.DictReader(open(A + "igfbp5_cells.csv"))}
def smooth_curve(d, nb=150, k=5):
    sel = np.abs(d["r_for_jsd"] - 0.5) < 0.15; o = np.argsort(d["t"][sel]); t = d["t"][sel][o]; f1 = d["x"][sel][o]; f2 = d["y"][sel][o]
    edges = np.linspace(t.min(), t.max(), nb+1); idx = np.clip(np.searchsorted(edges, t)-1, 0, nb-1)
    cx = np.array([np.median(f1[idx == b]) for b in range(nb) if np.any(idx == b)]); cy = np.array([np.median(f2[idx == b]) for b in range(nb) if np.any(idx == b)])
    ker = np.ones(k)/k; pad = k//2
    return (np.convolve(np.r_[cx[:pad][::-1], cx, cx[-pad:][::-1]], ker, "valid"), np.convolve(np.r_[cy[:pad][::-1], cy, cy[-pad:][::-1]], ker, "valid"))
def kde(vals, w=None, grid=np.linspace(0, 1, 200)):
    vals = np.asarray(vals); w = np.ones_like(vals) if w is None else np.asarray(w, float); w = w/w.sum()
    sd = np.sqrt(np.sum(w*(vals-np.sum(w*vals))**2)); n_eff = 1/np.sum(w**2); bw = max(0.9*sd*n_eff**(-0.2), 0.015)
    return grid, (w[None, :]*np.exp(-0.5*((grid[:, None]-vals[None, :])/bw)**2)).sum(1)/(bw*np.sqrt(2*np.pi))

fig = plt.figure(figsize=(6.5, 6.4))
def T(x, y, s, size=6.2, color=GREY, ha="center", va="center", **kw): return fig.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, **kw)
def plabel(x, y, s): fig.text(x, y, s, fontsize=11, fontweight="bold", va="top", ha="left")
def strip(rect):
    ax = fig.add_axes(rect); ax.set_xlim(-0.02, 1.02)
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    ax.set_yticks([]); ax.set_xticks([0, 0.5, 1]); ax.tick_params(axis="x", labelsize=5.4, length=2, pad=1); return ax
def chevron(x0, x1, y0, y1, color, label):
    n = 0.012; pts = [(x0, y0), (x1-n, y0), (x1, (y0+y1)/2), (x1-n, y1), (x0, y1), (x0+n, (y0+y1)/2)]
    fig.add_artist(Polygon(pts, closed=True, transform=fig.transFigure, fc="white", ec=color, lw=1.1)); T((x0+x1)/2, (y0+y1)/2, label, size=6.4, color=DARK)

# ================================================================ A  shared preprocessing (top half)
plabel(0.01, 0.99, "A"); T(0.42, 0.975, "Tissue-intrinsic coordinates", size=8, color=NAVY, fontweight="bold")
XS = (0.08, 0.31, 0.54); W = 0.22; H = 0.185
heads = ["Input sections", "Anatomical mapping", "Coordinates $(\\ell, r)$"]
for x0, lab in zip(XS, heads): T(x0 + W/2, 0.943, lab, size=6.6, color=DARK, fontweight="bold")
for x0 in XS[:-1]:
    fig.add_artist(FancyArrowPatch((x0 + W - 0.012, 0.943), (x0 + W + 0.012, 0.943), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=7, lw=0.8, color=GREY))
rows = [("D0", "Day 0\n(healthy)", DARK, 0.735), ("D9", "Day 9\n(DSS colitis)", DARK, 0.545)]
AX3 = {}
for key, lab, col, yb in rows:
    T(0.045, yb+0.09, lab, size=6.4, color=col, fontweight="bold")
    axs = [fig.add_axes((x, yb, W, H)) for x in XS]
    for L in ("ME", "SM", "MU"):
        d = SEC[key][L]
        axs[0].scatter(d["x"], d["y"], s=0.12, c=LMUTED[L], lw=0, alpha=0.9, rasterized=True)
        axs[1].scatter(d["x"], d["y"], s=0.12, c=LC[L], lw=0, alpha=0.75, rasterized=True)
        cx, cy = smooth_curve(d); axs[1].plot(cx, cy, color=CURVE, lw=0.4, path_effects=[pe.Stroke(linewidth=1.1, foreground="white"), pe.Normal()])
        axs[2].scatter(d["x"], d["y"], s=0.12, c=np.clip(d["r_for_jsd"], 0, 1), cmap=LCM[L], vmin=-0.18, vmax=1.0, lw=0, alpha=0.95, rasterized=True)
    for a in axs: a.set_aspect("equal"); a.axis("off")
    AX3[key] = axs[2]
# magnified mucosal inset: r runs from the lumen (0) to the outer boundary (1) within one layer
mu = SEC["D0"]["MU"]; tc = float(os.environ.get("TC", "0.85")); sel = np.abs(mu["t"] - tc) < 0.03
cx0, cy0 = np.median(mu["x"][sel]), np.median(mu["y"][sel]); hw = 200
INS = fig.add_axes((0.785, 0.745, 0.175, 0.175)); INS.set_xlim(cx0-hw, cx0+hw); INS.set_ylim(cy0-hw, cy0+hw); INS.set_aspect("equal")
for L in ("ME", "SM", "MU"):
    d = SEC["D0"][L]; m = (np.abs(d["x"]-cx0) < hw) & (np.abs(d["y"]-cy0) < hw)
    INS.scatter(d["x"][m], d["y"][m], s=7, c=np.clip(d["r_for_jsd"][m], 0, 1), cmap=LCM[L], vmin=-0.18, vmax=1.0, lw=0.15, edgecolors="white", rasterized=True)
# arrow direction = normal to the local mucosal medial curve (the direction along which r is measured), oriented toward increasing r
m = (np.abs(mu["x"]-cx0) < hw) & (np.abs(mu["y"]-cy0) < hw); rr = mu["r_for_jsd"][m]; XY = np.c_[mu["x"][m]-cx0, mu["y"][m]-cy0]
B = XY[rr >= np.quantile(rr, 0.85)]; tang = np.linalg.svd(B - B.mean(0))[2][0]; u = np.array([-tang[1], tang[0]])  # normal to the outer mucosal boundary in the window
grad = np.linalg.lstsq(np.c_[np.ones(len(rr)), XY], rr, rcond=None)[0][1:]
if grad @ u < 0: u = -u
proj = XY @ u; q10, q90 = np.quantile(rr, [0.03, 0.97]); s0 = np.median(proj[rr <= q10]); s1 = np.median(proj[rr >= q90]); lim = hw - 45
s0, s1 = max(s0, -lim), min(s1, lim); a0 = np.array([cx0, cy0]) + s0*u; a1 = np.array([cx0, cy0]) + s1*u
INS.annotate("", a1, a0, arrowprops=dict(arrowstyle="-|>", color=DARK, lw=1.3, mutation_scale=9, shrinkA=0, shrinkB=0), zorder=5)
v = np.array([-u[1], u[0]])
for p, lab, off in ((a0, "0", -1), (a1, "1", 1)):
    q = p + 24*off*u; q = np.clip(q, [cx0-hw+18, cy0-hw+18], [cx0+hw-18, cy0+hw-18])
    INS.text(q[0], q[1], lab, fontsize=6.2, ha="center", va="center", color=DARK, fontweight="bold", zorder=6, path_effects=[pe.Stroke(linewidth=1.8, foreground="white"), pe.Normal()])
qm = (a0 + a1)/2 + 26*v
INS.text(qm[0], qm[1], "r", fontsize=6.2, ha="center", va="center", color=DARK, fontweight="bold", style="italic", zorder=6, path_effects=[pe.Stroke(linewidth=1.8, foreground="white"), pe.Normal()])
INS.set_xticks([]); INS.set_yticks([])
for sp in INS.spines.values(): sp.set_color("#888888"); sp.set_linewidth(0.6)
T(0.785 + 0.0875, 0.927, "mucosa, magnified", size=5.6, color=DARK, va="bottom")
T(0.785 + 0.0875, 0.738, "r is computed within each layer\nlumen 0 \u2192 outer boundary 1", size=5.0, va="top")
# box on the Day 0 coordinate panel + connector to the inset
AX3["D0"].add_patch(Rectangle((cx0-hw, cy0-hw), 2*hw, 2*hw, fc="none", ec=DARK, lw=0.6, zorder=5))
bx, by = fig.transFigure.inverted().transform(AX3["D0"].transData.transform((cx0+hw, cy0)))
fig.add_artist(Line2D([bx, 0.785], [by, 0.745 + 0.0875], transform=fig.transFigure, color="#888888", lw=0.6, ls=(0, (2, 1.5))))
# legends under A (compact, below the inset)
LX = 0.80
T(LX, 0.700, "Layers", size=5.8, color=DARK, ha="left", fontweight="bold")
for i, L in enumerate(("MU", "SM", "ME")):
    fig.add_artist(Rectangle((LX, 0.676-i*0.020), 0.011, 0.010, transform=fig.transFigure, fc=LC[L], ec="none")); T(LX+0.016, 0.681-i*0.020, LN[L], size=5.4, ha="left")
fig.add_artist(Line2D([LX, LX+0.02], [0.618, 0.618], transform=fig.transFigure, color=CURVE, lw=0.9)); T(LX+0.025, 0.618, "medial curve", size=5.2, ha="left")
T(LX, 0.597, "Within-layer position r", size=5.8, color=DARK, ha="left", fontweight="bold")
for i, L in enumerate(("MU", "SM", "ME")):
    cax = fig.add_axes((LX, 0.567-i*0.019, 0.11, 0.012)); cax.imshow(np.linspace(0.15, 1, 100)[None, :], aspect="auto", cmap=LCM[L], extent=(0, 1, 0, 1)); cax.set_xticks([]); cax.set_yticks([])
T(LX, 0.524, "0  inner", size=5.2, ha="left", va="top"); T(LX+0.11, 0.524, "outer  1", size=5.2, ha="right", va="top")
# branch connector: coordinates feed both analyses
fig.add_artist(Line2D([0.515, 0.515], [0.535, 0.485], transform=fig.transFigure, color=GREY, lw=0.9))
fig.add_artist(Line2D([0.27, 0.76], [0.485, 0.485], transform=fig.transFigure, color=GREY, lw=0.9))
for x in (0.27, 0.76): fig.add_artist(FancyArrowPatch((x, 0.485), (x, 0.462), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=8, lw=0.9, color=GREY))

# ================================================================ B  population correspondence (bottom left)
plabel(0.01, 0.46, "B"); T(0.27, 0.448, "Population correspondence", size=8, color=NAVY, fontweight="bold")
T(0.27, 0.428, "anatomical similarity", size=5.6)
mu0 = SEC["D0"]["MU"]; mu9 = SEC["D9"]["MU"]
rA = mu0["r_for_jsd"][mu0["Tier3"] == "Fibro 2"]; r4 = mu9["r_for_jsd"][mu9["Tier3"] == "IAF 4"]; r2 = mu9["r_for_jsd"][mu9["Tier3"] == "IAF 2"]
g, yA = kde(rA); _, y4 = kde(r4); _, y2 = kde(r2); ymax = max(yA.max(), y4.max(), y2.max())*1.5
for k, (yc, col, nm, verdict) in enumerate(((y4, PURPLE, "IAF 4", "similar positions"), (y2, ORANGE, "IAF 2", "different positions"))):
    ax = strip((0.06 + k*0.235, 0.265, 0.2, 0.145)); ax.set_ylim(0, ymax)
    ax.fill_between(g, 0, yA, color=DARK, alpha=0.15, lw=0); ax.plot(g, yA, color=DARK, lw=1.0); ax.plot(g, yc, color=col, lw=1.2)
    ax.text(0.0, ymax*0.99, "Fibro 2 (Day 0)", fontsize=5.6, va="top", color=DARK, fontweight="bold"); ax.text(1.0, ymax*0.99, f"{nm} (Day 9)", fontsize=5.6, va="top", ha="right", color=col, fontweight="bold")
    ax.text(0.5, ymax*0.80, verdict, fontsize=5.2, va="top", ha="center", color=col)
    ax.set_xlabel("r within the mucosa", fontsize=5.4, labelpad=0)
# ranked dot chart: all Day 9 fibroblast-lineage candidates, ordered by pooled similarity (no numbers on the axis)
import csv as _csv
rows_ = list(_csv.reader(open(A + "pooled_sim_sub.csv"))); hdr = rows_[0][1:]; f2 = [r for r in rows_[1:] if r[0] == "Fibro 2"][0]
cands = sorted(((c, float(v)) for c, v in zip(hdr, f2[1:])), key=lambda t: -t[1])
ax = fig.add_axes((0.14, 0.07, 0.34, 0.14)); ax.set_xlim(0, 0.55); ax.set_ylim(len(cands)-0.4, -0.9)
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.set_yticks([]); ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5]); ax.tick_params(axis="x", labelsize=5.2, length=2, pad=1); ax.spines["bottom"].set_color("#999999")
ax.set_xlabel("anatomical similarity score to Fibro 2 (Day 0)", fontsize=5.6, labelpad=1)
hl = {"IAF 4": PURPLE, "IAF 2": ORANGE}
for i, (c, v) in enumerate(cands):
    col = hl.get(c, "#9a9a9a"); ax.plot([0, v], [i, i], color=col, lw=0.8, alpha=0.6); ax.plot(v, i, "o", ms=4.2 if c in hl else 3.2, color=col)
    ax.text(-0.015, i, c, fontsize=5.6 if c in hl else 5.2, ha="right", va="center", color=col, fontweight="bold" if c in hl else "normal")
ax.text(0, -0.9, "All Day 9 candidate populations, ranked", fontsize=5.8, ha="left", va="bottom", color=DARK)

# ================================================================ C  gene repositioning (bottom right)
plabel(0.5, 0.46, "C"); T(0.76, 0.448, "Gene repositioning", size=8, color=NAVY, fontweight="bold")
DELTAS = {}
T(0.76, 0.428, "transcript position relative to its cell population", size=5.6)
for k, (key, slide, ct, lab, col) in enumerate((("D0", D0[1], "Fibro 2", "Day 0 · Fibro 2", DARK), ("D9", D9[1], "IAF 4", "Day 9 · IAF 4", DARK))):
    d = SEC[key]["MU"]; m = d["Tier3"] == ct; r = d["r_for_jsd"][m]; cid = d["Cell_ID"][m]
    cnt = np.array([igf.get((slide, c), 0) for c in cid]); ex = cnt > 0; r_ex = r[ex]; w = cnt[ex]
    mean_cells = r.mean(); mean_tx = np.sum(w*r_ex)/w.sum()
    g, yc = kde(r); _, yt = kde(r_ex, w); ymx = max(yc.max(), yt.max())*1.5
    ax = strip((0.555 + k*0.225, 0.265, 0.19, 0.145)); ax.set_ylim(0, ymx)
    ax.fill_between(g, 0, yc, color="#aaaaaa", alpha=0.18, lw=0); ax.plot(g, yc, color="#777777", lw=1.0, ls="--"); ax.plot(g, yt, color=TEAL, lw=1.2)
    ax.plot([mean_cells]*2, [0, ymx*0.70], color="#777777", lw=0.8, ls=":"); ax.plot([mean_tx]*2, [0, ymx*0.70], color=TEAL, lw=0.8, ls=":")
    ax.annotate("", (mean_tx, ymx*0.75), (mean_cells, ymx*0.75), arrowprops=dict(arrowstyle="-|>", color=TEAL, lw=1.0, mutation_scale=7))
    ax.text((mean_cells+mean_tx)/2, ymx*0.79, "δ", fontsize=7, ha="center", va="bottom", color=TEAL, fontweight="bold", style="italic")
    ax.text(0.0, ymx*0.99, lab, fontsize=5.6, va="top", color=col, fontweight="bold")
    ax.set_xlabel("r within the mucosa", fontsize=5.4, labelpad=0)
    DELTAS[key] = mean_tx - mean_cells; print(key, ct, round(mean_tx-mean_cells, 3))
fig.add_artist(Line2D([0.565, 0.585], [0.222, 0.222], transform=fig.transFigure, color="#777777", lw=1.0, ls="--")); T(0.59, 0.222, "all cells of the population", size=5.2, ha="left")
fig.add_artist(Line2D([0.73, 0.75], [0.222, 0.222], transform=fig.transFigure, color=TEAL, lw=1.2)); T(0.755, 0.222, "transcript-weighted", size=5.2, ha="left")
T(0.76, 0.205, "dotted: the two mean positions;  δ = their difference", size=5.2)
T(0.676, 0.188, "Igfbp5", size=5.4, color="#555555", ha="right", style="italic", family="Arial"); T(0.680, 0.188, "\u00b7 Fibro 2 (Day 0) \u2192 IAF 4 (Day 9)", size=5.4, color="#555555", ha="left")
# 1-D contrast display: delta(Day 0) and delta(Day 9) on one axis
ax = fig.add_axes((0.575, 0.075, 0.38, 0.10)); ax.set_xlim(-0.07, 0.13); ax.set_ylim(0, 10); ax.axis("off")
ax.plot([-0.065, 0.125], [4.5, 4.5], color=DARK, lw=0.9); ax.plot([0, 0], [3.6, 5.4], color=DARK, lw=0.9); ax.text(0, 2.9, "0", fontsize=5.4, ha="center", va="top", color=DARK)
ax.text(-0.065, 2.9, "\u2190 toward lumen", fontsize=5, ha="left", va="top", color=GREY); ax.text(0.125, 2.9, "toward crypt base \u2192", fontsize=5, ha="right", va="top", color=GREY)
d0v, d9v = DELTAS["D0"], DELTAS["D9"]
ax.plot(d0v, 4.5, "o", ms=5.5, mfc="white", mec=DARK, mew=1.1, zorder=3); ax.text(d0v, 6.0, "\u03b4 (Day 0)", fontsize=5.6, ha="center", va="bottom", color=DARK)
ax.plot(d9v, 4.5, "o", ms=5.5, color=TEAL, zorder=3); ax.text(d9v, 6.0, "\u03b4 (Day 9)", fontsize=5.6, ha="center", va="bottom", color=TEAL)
ax.annotate("", (d9v-0.006, 4.5), (d0v+0.006, 4.5), arrowprops=dict(arrowstyle="-|>", color=NAVY, lw=1.4, mutation_scale=9), zorder=2)
ax.text((d0v+d9v)/2, 8.6, "Δ = δ(Day 9) − δ(Day 0)", fontsize=6.2, ha="center", va="bottom", color=NAVY, fontweight="bold")
ax.text((d0v+d9v)/2, 5.3, "Δ", fontsize=6.5, ha="center", va="bottom", color=NAVY, fontweight="bold", style="italic")

fig.savefig(OUT+"/fig1_framework_v11.pdf", dpi=300); fig.savefig(OUT+"/fig1_framework_v11.png", dpi=300); print("saved")
