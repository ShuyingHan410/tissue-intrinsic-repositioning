"""Abundance sensitivity, v2 (2026-11-01, SLAT replaces CAST). Definition as in methods_section.tex,
'Abundance sensitivity': within every origin row of every kept Day0 x Day9 pair, Spearman correlation between
the method's raw score for each candidate and that candidate's share of the fibroblast-lineage (Tier1 ==
Fibroblast) cells of the Day9 section; median over rows and pairs. Also: % of pairs in which IAF 2 or IAF 3 is
ranked first for each origin. Candidate sets / pair mapping identical to make_fig3_benchmark_v2.py.
Writes fig23_assets/abundance_sensitivity_by_method_v2.csv. Run: python3 abundance_sensitivity_v2.py"""
import os, re, csv, collections, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); OUT = os.path.dirname(os.path.abspath(__file__))
METHODS = ["morphogam", "moscot", "slat", "paste2", "stalign"]
d9a = sorted(f[:-4] for f in os.listdir(f"{ROOT}/all_slices") if "_D9_" in f and f.endswith(".csv"))
d9m = sorted(re.sub(r"^DSS9_all_per_slice_|_MU_xyt_mapped\.csv$", "", f) for f in os.listdir(f"{ROOT}/output_threelayer/per_slice_morphogam") if re.match(r"^DSS9_all_.*_MU_xyt_mapped\.csv$", f))
share = {}
for s in d9m:
    c = collections.Counter()
    with open(f"{ROOT}/all_slices/{s}.csv") as fh:
        rd = csv.reader(fh); hdr = next(rd); i1 = hdr.index("Tier1"); i2 = hdr.index("Tier2")
        for r in rd:
            if r[i1] == "Fibroblast": c[r[i2]] += 1
    tot = sum(c.values()); share[s] = {k: v / tot for k, v in c.items()}
def rankdata(x):
    x = np.asarray(x, float); o = np.argsort(x); r = np.empty(len(x)); r[o] = np.arange(1, len(x) + 1)
    for v in np.unique(x):
        m = x == v
        if m.sum() > 1: r[m] = r[m].mean()
    return r
def spearman(a, b):
    ra, rb = rankdata(a), rankdata(b)
    return np.nan if ra.std() == 0 or rb.std() == 0 else np.corrcoef(ra, rb)[0, 1]
def fl(v):
    try: return float(v)
    except ValueError: return np.nan
def read_mat(p):
    rows = list(csv.reader(open(p))); cols = rows[0][1:]
    return {r[0]: {c: fl(v) for c, v in zip(cols, r[1:])} for r in rows[1:]}
origins = ["Fibro 1", "Fibro 12", "Fibro 13", "Fibro 15", "Fibro 2", "Fibro 4", "Fibro 5", "Fibro 6", "Fibro 7", "IAF 2", "Pericyte 1", "Pericyte 2"]
out = []
for m in METHODS:
    pdir = f"{ROOT}/benchmark_out/day0_day9/{m}/partials"; rhos = []; rho6 = []; first = collections.defaultdict(lambda: [0, 0]); n = 0
    for f in sorted(os.listdir(pdir)):
        u = int(re.search(r"unit(\d+)", f).group(1))
        if m == "morphogam": d9 = d9m[u % 21]
        else:
            d9 = d9a[u % 25]
            if d9 not in d9m: continue
        M = read_mat(f"{pdir}/{f}"); n += 1
        for o, row in M.items():
            cands = [c for c in row if np.isfinite(row[c])]
            if len(cands) < 3: continue
            s = np.array([row[c] for c in cands]); a = np.array([share[d9].get(c, 0.0) for c in cands])
            r = spearman(s, a)
            if np.isfinite(r):
                rhos.append(r)
                if o == "Fibro 6": rho6.append(r)
            first[o][1] += 1
            if any(c in ("IAF 2", "IAF 3") for c in cands if np.isclose(row[c], s.max())): first[o][0] += 1
    rec = {"method": m, "pairs": n, "rows": len(rhos), "median_rho_all_rows": round(float(np.median(rhos)), 2), "median_rho_Fibro6": round(float(np.median(rho6)), 2)}
    for o in origins: rec[f"IAF23_first_{o}"] = round(100 * first[o][0] / first[o][1]) if first[o][1] else ""
    out.append(rec); print(rec["method"], rec["median_rho_all_rows"], rec["median_rho_Fibro6"], {o: rec[f"IAF23_first_{o}"] for o in ("Fibro 2", "Fibro 4", "Fibro 6", "Fibro 7")})
with open(f"{OUT}/fig23_assets/abundance_sensitivity_by_method_v2.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
print("wrote fig23_assets/abundance_sensitivity_by_method_v2.csv")
