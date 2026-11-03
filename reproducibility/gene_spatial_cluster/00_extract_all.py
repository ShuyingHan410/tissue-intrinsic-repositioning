#!/usr/bin/env python3
"""Step 0 (one job): pull ALL 990 columns (942 genes + 48 Blank probes) out of adata.h5ad
for the 31 slices that have MorphoGAM coordinates, and write:
  meta.csv.gz               one row per cell: Slice_ID, Cell_ID, Tier1, Tier2, day, mouse,
                            total_counts (REAL genes only) and blank_counts (sum over the
                            48 Blank barcodes - REQUIRED by the 2.813-sigma expression test)
  chunks/chunk_NN.csv.gz    Slice_ID, Cell_ID + CHUNK_SIZE gene columns
  genes_index.csv           gene -> chunk map (blanks flagged)
Reads the h5ad in 50k-row slabs, so memory stays ~2 GB regardless of file size."""
import h5py, numpy as np, pandas as pd, os, sys

ADATA  = os.environ.get("ADATA", "../../adata.h5ad")
COORD  = os.environ.get("COORD_DIR", "../../output_threelayer/per_slice_morphogam")
OUT    = os.environ.get("OUT", ".")
CHUNK  = int(os.environ.get("CHUNK_SIZE", 30))

slices = sorted({f.split("_per_slice_")[1].rsplit("_MU_xyt_mapped.csv", 1)[0]
                 for f in os.listdir(COORD) if f.endswith("_MU_xyt_mapped.csv")})
print(f"{len(slices)} coordinate slices", flush=True)

f   = h5py.File(ADATA, "r")
dec = lambda a: [x.decode() if isinstance(x, bytes) else x for x in a]
var = dec(f["raw/var/_index"][:])
real = [i for i, g in enumerate(var) if not g.startswith("Blank")]

sl_cat, sl_code = dec(f["obs/Slice_ID/categories"][:]), f["obs/Slice_ID/codes"][:]
mask = np.isin(sl_code, [i for i, s in enumerate(sl_cat) if s in slices])
print(f"cells kept: {mask.sum():,} / {len(mask):,}", flush=True)

t1c, t1 = dec(f["obs/Tier1/categories"][:]), f["obs/Tier1/codes"][:]
t2c, t2 = dec(f["obs/Tier2/categories"][:]), f["obs/Tier2/codes"][:]

blanks = [i for i, g in enumerate(var) if g.startswith("Blank")]
X = f["raw/X"]; n = X.shape[0]; STEP = 50_000
idxs, tots, blks, mats = [], [], [], []
for a in range(0, n, STEP):
    b = min(a + STEP, n); m = mask[a:b]
    if not m.any(): continue
    blk = X[a:b, :][m]
    idxs.append(np.arange(a, b)[m])
    tots.append(blk[:, real].sum(1))          # real transcripts only
    blks.append(blk[:, blanks].sum(1))        # the 48 Blank barcodes
    mats.append(blk)
    print(f"  {b:,}", flush=True)
idx = np.concatenate(idxs); M = np.vstack(mats).astype(np.int32)

sid = np.array([sl_cat[c] for c in sl_code[idx]])
meta = pd.DataFrame({"Slice_ID": sid, "Cell_ID": dec(f["obs/Cell_ID"][:][idx]),
                     "Tier1": [t1c[c] for c in t1[idx]], "Tier2": [t2c[c] for c in t2[idx]],
                     "day": np.where(pd.Series(sid).str.contains("_D9_"), "Day 9", "Day 0"),
                     "mouse": pd.Series(sid).str.replace(r"_slice_.*", "", regex=True),
                     "total_counts": np.concatenate(tots),
                     "blank_counts": np.concatenate(blks)})
meta.to_csv(os.path.join(OUT, "meta.csv.gz"), index=False)
print(f"meta.csv.gz written: {meta.shape}", flush=True)
assert "blank_counts" in meta.columns, "blank_counts missing - the expression test needs it"
print(f"  n_blank_barcodes = {len(blanks)}  (inclusion_criteria.R assumes this for the 2.813-sigma test)", flush=True)

os.makedirs(os.path.join(OUT, "chunks"), exist_ok=True)
gmap = []
for ci, s in enumerate(range(0, len(var), CHUNK), start=1):
    cols = list(range(s, min(s + CHUNK, len(var))))
    df = pd.DataFrame({"Slice_ID": sid, "Cell_ID": meta.Cell_ID})
    for c in cols: df[var[c]] = M[:, c]
    df.to_csv(os.path.join(OUT, "chunks", f"chunk_{ci:02d}.csv.gz"), index=False)
    gmap += [{"gene": var[c], "chunk": ci, "is_blank": var[c].startswith("Blank")} for c in cols]
    print(f"chunk_{ci:02d}: {len(cols)} genes", flush=True)
pd.DataFrame(gmap).to_csv(os.path.join(OUT, "genes_index.csv"), index=False)
print(f"DONE: {ci} chunks", flush=True)
