"""
Shared data-prep + aggregation for the Day0 half-split benchmark.

Positive-control design
-----------------------
Take ONE Day0 (Healthy) slice, split it into two stratified halves (A, B) within
each Tier2 cell type, and treat the halves as pseudo-replicates. Each benchmark
method (moscot / PASTE2 / STalign) produces a cell-by-cell correspondence between the
two halves; we aggregate it into a Tier2 x Tier2 matrix (rows = half A, cols = half B)
and row-normalize. Because the two halves are statistically identical, a correct
method should put the largest mass on the diagonal (each type -> itself), mirroring
the self-match rank-1 behaviour of the three-step pipeline.

Used by: benchmark_moscot_day0split.Rmd,
         benchmark_paste2_day0split.Rmd
"""
import os
import numpy as np
import pandas as pd

_PROXY_KEYS = (
    ("http_proxy", "HTTP_PROXY"),
    ("https_proxy", "HTTPS_PROXY"),
    ("all_proxy", "ALL_PROXY"),
    ("no_proxy", "NO_PROXY"),
)


def configure_proxy(params=None):
    """Apply HTTP(S) proxy env vars before network-using imports or downloads.

    Reads http_proxy, https_proxy, all_proxy, and no_proxy from *params* (or their
    UPPERCASE forms). Non-empty param values override the environment for that key.
    When only one case variant is set, both lower- and upper-case names are synced
    so pip, git, requests, and curl all see the same settings.
    """
    p = dict(params or {})
    applied = []
    for lower, upper in _PROXY_KEYS:
        val = p.get(lower)
        if val is None:
            val = p.get(upper)
        if val is not None and str(val).strip():
            val = str(val).strip()
        else:
            val = os.environ.get(upper) or os.environ.get(lower)
        if val:
            os.environ[lower] = val
            os.environ[upper] = val
            if lower in p or upper in p:
                applied.append(f"{lower}={val}")
    if applied:
        print("proxy configured:", ", ".join(applied))
    return applied

# NOTE: scanpy / h5py / scipy are imported lazily inside the functions that need
# them, so this module imports cleanly in a minimal env (e.g. the STalign env,
# which only needs the h5py-based obs loader and has no scanpy).

# ---------------------------------------------------------------------------
# Canonical cell-type sets, shared by ALL Day0 half-split benchmarks so every
# method is scored on the SAME types as benchmark_lineage_clean.Rmd.
#
# These are the INTERSECTION of the types present in BOTH data sources:
#   * the MorphoGAM three-layer output (output_threelayer/per_slice_morphogam)
#     that feeds the three-step method, and
#   * the raw adata.h5ad that moscot / PASTE2 / STalign read.
# Two filters shrink the type set the three-step method can show, so we match it:
#   (1) MorphoGAM ">30 cells per (slice x anatomical layer)" drops rare types
#       entirely (M cells, Epithelial (Clu+), IAE 3, FRC, IAF 1/3/4/5); then
#   (2) three_step()'s own min_cells>=30 per split-half drops the still-sparse
#       EEC + "Repair associated  (Arg1+)" (epithelial) and IAF 2 (fibroblast).
# Two benchmark DESIGNS keep different type sets, because min_cells>=30 counts over
# a different unit in each (a type is kept if it clears >=30 in a layer on BOTH sides):
#   * Benchmark A "halfsplit" - both sides are HALVES of one slice -> 7 epi / 11 fibro
#     (the set benchmark_lineage_clean.Rmd displays).
#   * Benchmark B "slicepair" - both sides are FULL slices (~2x cells) -> more types
#     clear the bar: +EEC +"Repair associated  (Arg1+)" (epi), +IAF 2 (fibro) = 9 / 12.
# Restricting the OT methods to the SAME per-design set is the like-for-like
# comparison. Order = crypt-axis base -> tip (matches the three-step heatmap axis).
LINEAGE_TYPES = {          # Benchmark A (half-split)
    "Fibroblast": ["Fibro 6", "Fibro 13", "Fibro 2", "Fibro 7", "Fibro 5", "Fibro 12",
                   "Fibro 1", "Fibro 15", "Fibro 4", "Pericyte 1", "Pericyte 2"],
    "Epithelial": ["Stem cells", "TA", "Colonocytes", "Goblet 1", "Goblet 2",
                   "IAE 1", "IAE 2"],
}
LINEAGE_TYPES_SLICEPAIR = {   # Benchmark B (two full slices)
    "Fibroblast": ["Fibro 6", "Fibro 13", "Fibro 2", "Fibro 7", "Fibro 5", "Fibro 12",
                   "Fibro 1", "Fibro 15", "Fibro 4", "Pericyte 1", "Pericyte 2", "IAF 2"],
    "Epithelial": ["Stem cells", "TA", "Colonocytes", "Goblet 1", "Goblet 2", "EEC",
                   "IAE 1", "IAE 2", "Repair associated  (Arg1+)"],
}

_TYPES_BY_DESIGN = {"halfsplit": LINEAGE_TYPES, "slicepair": LINEAGE_TYPES_SLICEPAIR}

# Day0->Day9 cross-condition (Fibroblast only) does NOT impose a fixed type set:
# each slice pair keeps whatever fibro/IAF types clear >=`day0day9_min_cells` cells on
# its own side (Day0 rows, Day9 cols), mirroring the MorphoGAM linkage pipeline where
# the candidate set varies per pair. This ordered list is used ONLY to sort output
# columns for display; it never drops a present type.
DAY0DAY9_ORDER = ["Fibro 6", "Fibro 13", "Fibro 2", "IAF 2", "IAF 3", "IAF 4",
                  "Fibro 7", "Fibro 4", "Fibro 5", "Fibro 12", "Fibro 1", "Fibro 15",
                  "Pericyte 1", "Pericyte 2", "IAF 1", "IAF 5", "FRC"]
DAY0DAY9_ORIGINS = ["Fibro 2", "Fibro 6", "Fibro 7"]   # the 3 Day0 ground-truth rows

# Ground-truth lineage trajectory (ordinal positions along the crypt axis) used by
# the Concordance Index. Only these on-axis types carry a distance relationship;
# every other type is "off-axis" (identity constraint only). Same for both designs.
#   Epithelial crypt: Stem -> TA -> Colonocytes
#   Fibroblast (Cadinu 2024): Fibro 6 (base) -> Fibro 13 (mid) -> Fibro 2 (top)
LINEAGE_AXIS = {
    "Epithelial": {"Stem cells": 0, "TA": 1, "Colonocytes": 2},
    "Fibroblast": {"Fibro 6": 0, "Fibro 13": 1, "Fibro 2": 2},
}


def type_set(params):
    """Ordered type list for this run's design + lineage (None if lineage == 'all').

    design defaults to 'halfsplit' (Benchmark A) when not given.
    For 'day0day9' returns None: no fixed input restriction, each pair keeps its own
    types (see load_cross_pair, which applies a per-side min-cell filter instead).
    """
    if str(params.get("design", "halfsplit")) == "day0day9":
        return None
    table = _TYPES_BY_DESIGN.get(str(params.get("design", "halfsplit")), LINEAGE_TYPES)
    return table.get(str(params["lineage"]))


def lineage_axis(params):
    """Ordinal-position dict for the run's lineage ({} if none defined)."""
    return LINEAGE_AXIS.get(str(params.get("lineage", "")), {})


def _todense(x):
    import scipy.sparse as sp
    return x.toarray() if sp.issparse(x) else np.asarray(x)


def _stratified_half_subsample(labels, n_cells, seed, nmax):
    """Return boolean-ish 'half' array (A/B) and capped index arrays for each half."""
    rng = np.random.default_rng(int(seed))
    half = np.empty(n_cells, dtype=object)
    for t in np.unique(labels):
        idx = np.where(labels == t)[0]
        rng.shuffle(idx)
        h = len(idx) // 2
        half[idx[:h]] = "A"
        half[idx[h:]] = "B"
    def cap(mask):
        idx = np.where(mask)[0]
        if len(idx) > nmax:
            idx = np.sort(rng.choice(idx, nmax, replace=False))
        return idx
    return half, cap(half == "A"), cap(half == "B")


def _read_obs_col(obs, name):
    """Read one .obs column from an open h5ad (h5py), decoding categoricals."""
    import h5py
    node = obs[name]
    if isinstance(node, h5py.Group):                 # categorical
        cats = node["categories"][:]
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in cats])
        codes = node["codes"][:]
        return cats[codes.clip(min=0)]
    vals = node[:]
    if vals.dtype.kind == "S":
        vals = np.array([v.decode() for v in vals])
    return vals


def load_half_split(params):
    """Load one Healthy slice, stratified-split into halves A/B, subsample.

    params (an R list passed via reticulate -> dict) must contain:
      adata_path, sample_type, slice_id ("FIRST" to auto-pick), lineage
      ("all" or a Tier1 value), seed, max_cells_per_half
    Returns: (A, B, types, slice_id)
      A, B   : AnnData halves (with .obsm['spatial'] = [x, y])
      types  : sorted list of Tier2 labels present (for a stable heatmap order)
      slice_id : the slice actually used
    """
    import scanpy as sc
    p = dict(params)
    adata = sc.read_h5ad(p["adata_path"], backed="r")
    obs = adata.obs
    sids = obs["Slice_ID"].astype(str)

    m = (obs["Sample_type"].astype(str) == str(p["sample_type"])).values
    slice_id = str(p["slice_id"])
    if slice_id == "FIRST":                       # busiest slice of this condition
        slice_id = sids[m].value_counts().index[0]
    m = m & (sids == slice_id).values

    lineage = str(p["lineage"])
    if lineage.lower() != "all":                  # restrict to one Tier1 lineage
        m = m & (obs["Tier1"].astype(str) == lineage).values
        keep = type_set(p)                        # input-level fairness: keep only the
        if keep is not None:                      # types the three-step method computes on,
            m = m & obs["Tier2"].astype(str).isin(keep).values  # so no OT mass leaks to
                                                  # excluded types (per-design: A=7/11, B=9/12)

    sub = adata[m].to_memory()                    # loads ONLY the selected cells
    if sub.raw is not None:
        sub.raw = None                            # avoid materializing the full .raw
    sub.obs["Tier2"] = sub.obs["Tier2"].astype(str)
    sub.obsm["spatial"] = sub.obs[["x", "y"]].to_numpy(dtype="float32")

    half, ia, ib = _stratified_half_subsample(
        sub.obs["Tier2"].values, sub.n_obs, p["seed"], int(p["max_cells_per_half"]))
    sub.obs["half"] = half
    A = sub[ia].copy()
    B = sub[ib].copy()
    types = sorted(set(sub.obs["Tier2"]))
    print(f"slice={slice_id}  lineage={lineage}  "
          f"A={A.n_obs} cells  B={B.n_obs} cells  types={len(types)}")
    return A, B, types, slice_id


def load_half_split_obs(params):
    """h5py-only variant: returns (A_df, B_df, types, slice_id) with x/y/Tier2.

    For position-only methods (STalign) that need no expression and run in an env
    without scanpy/anndata. A_df/B_df are pandas DataFrames with columns x, y, Tier2.
    """
    import h5py
    p = dict(params)
    f = h5py.File(p["adata_path"], "r")
    obs = f["obs"]
    df = pd.DataFrame({
        "Sample_type": _read_obs_col(obs, "Sample_type").astype(str),
        "Slice_ID":    _read_obs_col(obs, "Slice_ID").astype(str),
        "Tier1":       _read_obs_col(obs, "Tier1").astype(str),
        "Tier2":       _read_obs_col(obs, "Tier2").astype(str),
        "x":           np.asarray(_read_obs_col(obs, "x"), dtype=float),
        "y":           np.asarray(_read_obs_col(obs, "y"), dtype=float),
    })
    f.close()

    m = df["Sample_type"] == str(p["sample_type"])
    slice_id = str(p["slice_id"])
    if slice_id == "FIRST":
        slice_id = df.loc[m, "Slice_ID"].value_counts().index[0]
    m = m & (df["Slice_ID"] == slice_id)
    lineage = str(p["lineage"])
    if lineage.lower() != "all":
        m = m & (df["Tier1"] == lineage)
        keep = type_set(p)                        # input-level fairness (see load_half_split)
        if keep is not None:
            m = m & df["Tier2"].isin(keep)
    sub = df[m].reset_index(drop=True)

    half, ia, ib = _stratified_half_subsample(
        sub["Tier2"].values, len(sub), p["seed"], int(p["max_cells_per_half"]))
    A = sub.iloc[ia].reset_index(drop=True)
    B = sub.iloc[ib].reset_index(drop=True)
    types = sorted(set(sub["Tier2"]))
    print(f"slice={slice_id}  lineage={lineage}  "
          f"A={len(A)} cells  B={len(B)} cells  types={len(types)}")
    return A, B, types, slice_id


def all_slice_ids(params):
    """All Slice_ID values for the chosen Sample_type (e.g. all 10 Healthy Day0 slices).

    Optionally restricted to slices that actually contain the chosen lineage. Uses
    h5py so it works in any env. Returns a sorted list of slice ids.
    """
    import h5py
    p = dict(params)
    f = h5py.File(p["adata_path"], "r")
    obs = f["obs"]
    st = _read_obs_col(obs, "Sample_type").astype(str)
    sid = _read_obs_col(obs, "Slice_ID").astype(str)
    t1 = _read_obs_col(obs, "Tier1").astype(str)
    f.close()
    m = st == str(p["sample_type"])
    lineage = str(p.get("lineage", "all"))
    if lineage.lower() != "all":
        m = m & (t1 == lineage)
    return sorted(set(sid[m]))


# --------------------------------------------------------------------------- #
# Benchmark B (slice-pair): two DIFFERENT full Day0 slices as A and B.
# --------------------------------------------------------------------------- #
def _cap(sub_n, seed, nmax):
    rng = np.random.default_rng(int(seed))
    if sub_n > nmax:
        return np.sort(rng.choice(sub_n, nmax, replace=False))
    return np.arange(sub_n)


def load_slice_pair(params):
    """Load slice_a as A and slice_b as B (both FULL slices, capped to
    max_cells_per_half), same lineage, input-filtered to the design's type set.
    Returns (A, B, types, tag) mirroring load_half_split."""
    import scanpy as sc
    p = dict(params)
    adata = sc.read_h5ad(p["adata_path"], backed="r")
    obs = adata.obs
    sids = obs["Slice_ID"].astype(str).values
    base = (obs["Sample_type"].astype(str) == str(p["sample_type"])).values
    lineage = str(p["lineage"])
    if lineage.lower() != "all":
        base = base & (obs["Tier1"].astype(str) == lineage).values
        keep = type_set(p)
        if keep is not None:
            base = base & obs["Tier2"].astype(str).isin(keep).values

    def one(sid):
        sub = adata[base & (sids == str(sid))].to_memory()
        if sub.raw is not None:
            sub.raw = None
        sub.obs["Tier2"] = sub.obs["Tier2"].astype(str)
        sub.obsm["spatial"] = sub.obs[["x", "y"]].to_numpy(dtype="float32")
        idx = _cap(sub.n_obs, p["seed"], int(p["max_cells_per_half"]))
        return sub[idx].copy()

    A = one(p["slice_a"]); B = one(p["slice_b"])
    types = sorted(set(A.obs["Tier2"]) | set(B.obs["Tier2"]))
    print(f"pair A={p['slice_a']}({A.n_obs}) B={p['slice_b']}({B.n_obs})  "
          f"lineage={lineage}  types={len(types)}", flush=True)
    return A, B, types, f"{p['slice_a']}__{p['slice_b']}"


def load_slice_pair_obs(params):
    """h5py-only slice-pair loader (STalign). Returns (A_df, B_df, types, tag)."""
    import h5py
    p = dict(params)
    f = h5py.File(p["adata_path"], "r")
    obs = f["obs"]
    df = pd.DataFrame({
        "Sample_type": _read_obs_col(obs, "Sample_type").astype(str),
        "Slice_ID":    _read_obs_col(obs, "Slice_ID").astype(str),
        "Tier1":       _read_obs_col(obs, "Tier1").astype(str),
        "Tier2":       _read_obs_col(obs, "Tier2").astype(str),
        "x":           np.asarray(_read_obs_col(obs, "x"), dtype=float),
        "y":           np.asarray(_read_obs_col(obs, "y"), dtype=float),
    })
    f.close()
    base = df["Sample_type"] == str(p["sample_type"])
    lineage = str(p["lineage"])
    if lineage.lower() != "all":
        base = base & (df["Tier1"] == lineage)
        keep = type_set(p)
        if keep is not None:
            base = base & df["Tier2"].isin(keep)

    def one(sid):
        sub = df[base & (df["Slice_ID"] == str(sid))].reset_index(drop=True)
        return sub.iloc[_cap(len(sub), p["seed"], int(p["max_cells_per_half"]))].reset_index(drop=True)

    A = one(p["slice_a"]); B = one(p["slice_b"])
    types = sorted(set(A["Tier2"]) | set(B["Tier2"]))
    print(f"pair A={p['slice_a']}({len(A)}) B={p['slice_b']}({len(B)})  "
          f"lineage={lineage}  types={len(types)}", flush=True)
    return A, B, types, f"{p['slice_a']}__{p['slice_b']}"


# --------------------------------------------------------------------------- #
# Day0 -> Day9 cross-condition: A = a `sample_type_a` (Healthy/Day0) slice,
# B = a `sample_type_b` (DSS9/Day9) slice. Rows = Day0 fibro types, cols = Day9
# fibro/IAF types. NO fixed type set: each side keeps its own Tier2 types with
# >= day0day9_min_cells cells, so the candidate set varies per pair (matches the
# MorphoGAM linkage pipeline). max_cells_per_half caps total cells per side.
# --------------------------------------------------------------------------- #
def _min_cell_keep(labels, min_cells):
    vc = pd.Series(labels).value_counts()
    return set(vc[vc >= int(min_cells)].index)


def load_cross_pair(params):
    import scanpy as sc
    p = dict(params); lineage = str(p["lineage"]); mc = int(p.get("day0day9_min_cells", 30))
    adata = sc.read_h5ad(p["adata_path"], backed="r"); obs = adata.obs
    sids = obs["Slice_ID"].astype(str).values; st = obs["Sample_type"].astype(str).values
    t1 = obs["Tier1"].astype(str).values
    def one(sid, sample_type):
        base = (st == str(sample_type)) & (sids == str(sid))
        if lineage.lower() != "all":
            base = base & (t1 == lineage)
        sub = adata[base].to_memory()
        if sub.raw is not None: sub.raw = None
        sub.obs["Tier2"] = sub.obs["Tier2"].astype(str)
        keep = _min_cell_keep(sub.obs["Tier2"].values, mc)          # per-pair type filter
        sub = sub[sub.obs["Tier2"].isin(keep)].copy()
        sub.obsm["spatial"] = sub.obs[["x", "y"]].to_numpy(dtype="float32")
        idx = _cap(sub.n_obs, p["seed"], int(p["max_cells_per_half"]))
        return sub[idx].copy()
    A = one(p["slice_a"], p.get("sample_type_a", "Healthy"))
    B = one(p["slice_b"], p.get("sample_type_b", "DSS9"))
    types = sorted(set(A.obs["Tier2"]) | set(B.obs["Tier2"]))
    print(f"cross A={p['slice_a']}({A.n_obs}) B={p['slice_b']}({B.n_obs}) "
          f"lineage={lineage} types={len(types)}", flush=True)
    return A, B, types, f"{p['slice_a']}__{p['slice_b']}"


def load_cross_pair_obs(params):
    """h5py-only cross-condition loader (STalign). Returns (A_df, B_df, types, tag)."""
    import h5py
    p = dict(params); lineage = str(p["lineage"]); mc = int(p.get("day0day9_min_cells", 30))
    f = h5py.File(p["adata_path"], "r"); obs = f["obs"]
    df = pd.DataFrame({
        "Sample_type": _read_obs_col(obs, "Sample_type").astype(str),
        "Slice_ID":    _read_obs_col(obs, "Slice_ID").astype(str),
        "Tier1":       _read_obs_col(obs, "Tier1").astype(str),
        "Tier2":       _read_obs_col(obs, "Tier2").astype(str),
        "x":           np.asarray(_read_obs_col(obs, "x"), dtype=float),
        "y":           np.asarray(_read_obs_col(obs, "y"), dtype=float)})
    f.close()
    def one(sid, sample_type):
        base = (df["Sample_type"] == str(sample_type)) & (df["Slice_ID"] == str(sid))
        if lineage.lower() != "all":
            base = base & (df["Tier1"] == lineage)
        sub = df[base].reset_index(drop=True)
        keep = _min_cell_keep(sub["Tier2"].values, mc)
        sub = sub[sub["Tier2"].isin(keep)].reset_index(drop=True)
        return sub.iloc[_cap(len(sub), p["seed"], int(p["max_cells_per_half"]))].reset_index(drop=True)
    A = one(p["slice_a"], p.get("sample_type_a", "Healthy"))
    B = one(p["slice_b"], p.get("sample_type_b", "DSS9"))
    types = sorted(set(A["Tier2"]) | set(B["Tier2"]))
    return A, B, types, f"{p['slice_a']}__{p['slice_b']}"


def load_pair(params):
    """Design dispatcher: half-split (A), slice-pair (B), or day0day9 cross-condition."""
    d = str(params.get("design", "halfsplit"))
    if d == "day0day9":  return load_cross_pair(params)
    if d == "slicepair": return load_slice_pair(params)
    return load_half_split(params)


def load_pair_obs(params):
    """h5py-only design dispatcher for STalign."""
    d = str(params.get("design", "halfsplit"))
    if d == "day0day9":  return load_cross_pair_obs(params)
    if d == "slicepair": return load_slice_pair_obs(params)
    return load_half_split_obs(params)


def all_units(params):
    """Work items to pool over: slice ids (half-split), all C(n,2) Day0 pairs
    (slice-pair), or the full Day0 x Day9 grid (day0day9)."""
    d = str(params.get("design", "halfsplit"))
    if d == "day0day9":
        a = all_slice_ids({**params, "sample_type": params.get("sample_type_a", "Healthy")})
        b = all_slice_ids({**params, "sample_type": params.get("sample_type_b", "DSS9")})
        return [(x, y) for x in a for y in b]
    slices = all_slice_ids(params)
    if d == "slicepair":
        return [(slices[i], slices[j])
                for i in range(len(slices)) for j in range(i + 1, len(slices))]
    return list(slices)


def make_typematrix(weight, types_A, types_B, display_types=None, row_normalize=True):
    """Build a (restricted, row-normalized) Tier2 x Tier2 DataFrame in memory."""
    M = couple_to_typematrix(weight, types_A, types_B)
    if display_types is not None:
        r = [t for t in display_types if t in M.index]
        c = [t for t in display_types if t in M.columns]
        M = M.loc[r, c]
    if row_normalize:
        rs = M.sum(axis=1).replace(0, np.nan)
        M = M.div(rs, axis=0)
    return M


def pool_matrices(mats, display_types=None):
    """Average a list of Tier2 x Tier2 matrices over slices (mirrors mean_matrix in R).

    Cells missing in a slice are skipped; each (row, col) is averaged over the slices
    where it is present.
    """
    mats = [m for m in mats if m is not None and m.size]
    if not mats:
        return pd.DataFrame()
    rows = sorted(set().union(*[m.index for m in mats]))
    cols = sorted(set().union(*[m.columns for m in mats]))
    if display_types is not None:
        rows = [t for t in display_types if t in rows]
        cols = [t for t in display_types if t in cols]
    acc = pd.DataFrame(0.0, index=rows, columns=cols)
    cnt = pd.DataFrame(0.0, index=rows, columns=cols)
    for m in mats:
        mm = m.reindex(index=rows, columns=cols)
        present = mm.notna()
        acc = acc.add(mm.fillna(0.0))
        cnt = cnt.add(present.astype(float))
    pooled = acc / cnt.replace(0, np.nan)
    return pooled


def couple_to_typematrix(weight, types_A, types_B):
    """Aggregate a cell-by-cell coupling/affinity (nA x nB) into Tier2 x Tier2.

    rows = half A types, cols = half B types.
    """
    dfA = pd.get_dummies(pd.Series(np.asarray(types_A)))   # nA x |TA|
    dfB = pd.get_dummies(pd.Series(np.asarray(types_B)))   # nB x |TB|
    M = dfA.T.values @ np.asarray(weight) @ dfB.values
    return pd.DataFrame(M, index=list(dfA.columns), columns=list(dfB.columns))


def save_typematrix(M, out_dir, name, display_types=None, row_normalize=True):
    """Optionally restrict to display_types, row-normalize, write CSV."""
    os.makedirs(out_dir, exist_ok=True)
    if display_types is not None:
        r = [t for t in display_types if t in M.index]
        c = [t for t in display_types if t in M.columns]
        M = M.loc[r, c]
    if row_normalize:
        rs = M.sum(axis=1).replace(0, np.nan)
        M = M.div(rs, axis=0)
    path = os.path.join(out_dir, name)
    M.to_csv(path)
    print("wrote", path, "shape", tuple(M.shape))
    return path
