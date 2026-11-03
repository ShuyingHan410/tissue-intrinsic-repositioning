#!/usr/bin/env python
"""Standalone driver for the Day0 half-split benchmark (no R / rmarkdown needed).

Runs one method (moscot | paste2 | stalign | cast) over BOTH lineages
(Fibroblast, Epithelial), pooling all 10 Day0 (Healthy) slices, and writes
<method>_typematrix_<lineage>.csv into --out. Uses bench_common, so the shared
7/11 input-cell filter (LINEAGE_TYPES) and slice pooling are applied identically
to every method. Plot the CSVs afterward with bench_heatmap.R for the heatmaps.

Examples (on FASRC):
  python run_bench.py --method paste2 --adata /path/adata.h5ad \
         --out /path/benchmark_out/paste2 --max_cells 4000
  python run_bench.py --method cast   --adata /path/adata.h5ad \
         --out /path/benchmark_out/cast   --max_cells 4000
"""
import argparse
import gc
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # find bench_common
import numpy as np
import pandas as pd
import bench_common as bc


# --------------------------------------------------------------------------- #
# Per-method, per-slice implementations. Each takes a params dict `p` (with
# slice_id + lineage set) and returns a Tier2 x Tier2 DataFrame (rows = half A,
# cols = half B), restricted to LINEAGE_TYPES[lineage] and row-normalized.
# --------------------------------------------------------------------------- #
def m_moscot(p):
    import scanpy as sc
    import anndata as adn
    from moscot.problems.space import AlignmentProblem
    A, B, _, _ = bc.load_pair(p)
    comb = adn.concat([A, B], label="batch", keys=["A", "B"], index_unique="-")
    comb.obsm["spatial"] = comb.obs[["x", "y"]].to_numpy(dtype="float32")
    if "X_pca" not in comb.obsm:
        sc.pp.pca(comb, n_comps=50)
    ap = (AlignmentProblem(adata=comb)
          .prepare(batch_key="batch", spatial_key="spatial", joint_attr="X_pca")
          .solve(alpha=0.5, epsilon=1e-2))
    (src, tgt), sol = list(ap.solutions.items())[0]
    pi = np.asarray(sol.transport_matrix)
    t_src = comb.obs["Tier2"][comb.obs["batch"] == src].astype(str).values
    t_tgt = comb.obs["Tier2"][comb.obs["batch"] == tgt].astype(str).values
    M = bc.make_typematrix(pi, t_src, t_tgt, display_types=bc.type_set(p),
                           row_normalize=p.get("row_normalize", False))
    return M if src == "A" else M.T


def m_paste2(p):
    try:
        from paste2 import PASTE2
        def align(A, B):
            return np.asarray(PASTE2.partial_pairwise_align(A, B, s=0.99))
    except Exception as e:
        print("PASTE2 unavailable, falling back to original PASTE:", e)
        import paste as pst
        def align(A, B):
            return np.asarray(pst.pairwise_align(A, B))
    A, B, _, _ = bc.load_pair(p)
    A.obsm["spatial"] = A.obs[["x", "y"]].to_numpy(dtype="float32")
    B.obsm["spatial"] = B.obs[["x", "y"]].to_numpy(dtype="float32")
    tA = A.obs["Tier2"].astype(str).values
    tB = B.obs["Tier2"].astype(str).values
    pi = align(A, B)
    M = bc.make_typematrix(pi, tA, tB, display_types=bc.type_set(p),
                           row_normalize=p.get("row_normalize", False))
    del A, B, pi
    gc.collect()
    return M


def m_stalign(p):
    from sklearn.neighbors import NearestNeighbors
    from STalign import STalign
    dx = float(p.get("raster_dx", 30))
    niter = int(p.get("niter", 200))
    device = "cpu"
    A, B, _, _ = bc.load_pair_obs(p)          # x/y/Tier2 only (no scanpy)
    A_xy = A[["x", "y"]].to_numpy(float)
    B_xy = B[["x", "y"]].to_numpy(float)
    xI, yI = A_xy[:, 0], A_xy[:, 1]
    xJ, yJ = B_xy[:, 0], B_xy[:, 1]
    XI, YI, I, _ = STalign.rasterize(xI, yI, dx=dx)
    XJ, YJ, J, _ = STalign.rasterize(xJ, yJ, dx=dx)
    I = I / I.mean(); J = J / J.mean()
    out = STalign.LDDMM([YI, XI], I, [YJ, XJ], J, niter=niter, device=device)
    Aaff, v, xv = out["A"], out["v"], out["xv"]
    ptsI = np.stack([yI, xI], axis=-1)
    tptsI = np.asarray(STalign.transform_points_source_to_target(xv, v, Aaff, ptsI))
    alignedA = np.stack([tptsI[:, 1], tptsI[:, 0]], axis=-1)
    nn_scale = np.median(NearestNeighbors(n_neighbors=1).fit(B_xy)
                         .kneighbors(alignedA)[0][:, 0])
    sigma = nn_scale * float(p.get("sigma_mult", 2))
    aa = (alignedA ** 2).sum(1)[:, None]; bb = (B_xy ** 2).sum(1)[None, :]
    D2 = np.maximum(aa + bb - 2.0 * alignedA @ B_xy.T, 0.0)
    K = np.exp(-D2 / (2.0 * sigma ** 2))
    u = np.ones(len(alignedA)); w = np.ones(len(B_xy))
    for _ in range(300):
        u = 1.0 / (K @ w + 1e-12); w = 1.0 / (K.T @ u + 1e-12)
    P = u[:, None] * K * w[None, :]
    return bc.make_typematrix(P, A["Tier2"].astype(str).values,
                              B["Tier2"].astype(str).values,
                              display_types=bc.type_set(p),
                              row_normalize=p.get("row_normalize", False))


def m_cast(p):
    import torch
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import normalize
    from CAST import CAST_MARK
    knn = int(p.get("knn_phys", 50))
    device = "cuda" if torch.cuda.is_available() else "cpu"
    A, B, _, _ = bc.load_pair(p)
    # Coordinates must stay on CPU: CAST_MARK builds the Delaunay graph with
    # scipy (calls .numpy()), which cannot convert a CUDA tensor. Only the
    # expression tensors go to the GPU for the GNN embedding.
    coords = {"A": torch.tensor(A.obsm["spatial"], dtype=torch.float32),
              "B": torch.tensor(B.obsm["spatial"], dtype=torch.float32)}
    expr = {"A": torch.tensor(bc._todense(A.X), dtype=torch.float32).to(device),
            "B": torch.tensor(bc._todense(B.X), dtype=torch.float32).to(device)}
    # Unique embedding dir per work unit: in slicepair mode slice_id is always
    # "FIRST", so concurrent array tasks must not share cast_mark/ or they race.
    tag = str(p.get("_urepr") or p.get("slice_id", "FIRST"))
    tag = "".join(c if (c.isalnum() or c in "._-") else "_" for c in tag)
    mark_dir = os.path.join(p["out_dir"], "cast_mark", tag)
    os.makedirs(mark_dir, exist_ok=True)
    embed = CAST_MARK(coords, expr, mark_dir)
    embA = normalize(embed["A"].detach().cpu().numpy())
    embB = normalize(embed["B"].detach().cpu().numpy())
    k = int(min(knn, A.n_obs))
    nn = NearestNeighbors(n_neighbors=k).fit(A.obsm["spatial"])
    _, cand = nn.kneighbors(B.obsm["spatial"])
    match = np.array([cand[j][(embB[j] @ embA[cand[j]].T).argmax()]
                      for j in range(B.n_obs)])
    P = np.zeros((A.n_obs, B.n_obs), dtype=np.float32)
    P[match, np.arange(B.n_obs)] = 1.0
    return bc.make_typematrix(P, A.obs["Tier2"].astype(str).values,
                              B.obs["Tier2"].astype(str).values,
                              display_types=bc.type_set(p),
                              row_normalize=p.get("row_normalize", False))


def m_moscot_ub(p):
    """moscot fused GW with UNBALANCED marginals (tau_a = tau_b = 0.9): lets population mass be
    created/destroyed between conditions, the appropriate relaxation for cross-condition pairs
    where populations expand (IAFs) or vanish. Everything else identical to m_moscot."""
    import scanpy as sc
    import anndata as adn
    from moscot.problems.space import AlignmentProblem
    A, B, _, _ = bc.load_pair(p)
    comb = adn.concat([A, B], label="batch", keys=["A", "B"], index_unique="-")
    comb.obsm["spatial"] = comb.obs[["x", "y"]].to_numpy(dtype="float32")
    if "X_pca" not in comb.obsm:
        sc.pp.pca(comb, n_comps=50)
    tau = float(p.get("tau", 0.9))
    ap = (AlignmentProblem(adata=comb)
          .prepare(batch_key="batch", spatial_key="spatial", joint_attr="X_pca")
          .solve(alpha=0.5, epsilon=1e-2, tau_a=tau, tau_b=tau))
    (src, tgt), sol = list(ap.solutions.items())[0]
    pi = np.asarray(sol.transport_matrix)
    t_src = comb.obs["Tier2"][comb.obs["batch"] == src].astype(str).values
    t_tgt = comb.obs["Tier2"][comb.obs["batch"] == tgt].astype(str).values
    M = bc.make_typematrix(pi, t_src, t_tgt, display_types=bc.type_set(p),
                           row_normalize=p.get("row_normalize", False))
    return M if src == "A" else M.T


def m_slat(p):
    """SLAT (scSLAT, Xia et al. 2023): graph-adversarial alignment designed for heterogeneous /
    cross-condition section pairs. Pipeline: kNN spatial graph per section -> shared embedding by
    run_SLAT -> for every Day9 (B) cell, its top-k Day0 (A) neighbours in embedding space
    (spatial_match with smooth=False, so no assumption of a shared coordinate frame). The
    coupling is SOFT: each B cell spreads unit mass over its top-k A cells in proportion to the
    cosine similarity, which avoids the integer-count ties of a hard one-to-one match.
    Features: adata.X as stored (log-normalised; feature='raw' in load_anndatas) because
    bc.load_pair drops .raw, so SLAT's DPCA/HVG branches (which expect counts) are not usable."""
    from scSLAT.model import Cal_Spatial_Net, load_anndatas, run_SLAT, spatial_match
    k_graph = int(p.get("slat_k", 10)); top = int(p.get("slat_top", 20)); feat = str(p.get("slat_feature", "raw"))
    A, B, _, _ = bc.load_pair(p)
    for ad in (A, B):
        ad.obsm["spatial"] = ad.obs[["x", "y"]].to_numpy(dtype="float32")
        if hasattr(ad.X, "toarray"): ad.X = ad.X.toarray()
        Cal_Spatial_Net(ad, k_cutoff=k_graph, model="KNN")
    edges, features = load_anndatas([A, B], feature=feat, check_order=False)   # scSLAT only checks n_A >= n_B; irrelevant here
    embd0, embd1, _ = run_SLAT(features, edges)                     # embd0: A cells, embd1: B cells
    out = spatial_match([embd0, embd1], adatas=[A, B], reorder=False, smooth=False, smooth_range=top)
    order, sim = np.asarray(out[1]), np.asarray(out[2])             # (n_B, top): A indices and cosine similarities per B cell
    w = np.clip(sim, 0, None); w = w / np.maximum(w.sum(1, keepdims=True), 1e-12)
    P = np.zeros((A.n_obs, B.n_obs), dtype=np.float32)
    for j in range(B.n_obs):
        P[order[j], j] += w[j]
    return bc.make_typematrix(P, A.obs["Tier2"].astype(str).values,
                              B.obs["Tier2"].astype(str).values,
                              display_types=bc.type_set(p),
                              row_normalize=p.get("row_normalize", False))


def m_cast_stack(p):
    """CAST run as designed for cross-sample use (Tang et al. 2024): CAST_MARK embeddings, then
    CAST_STACK registration of the Day9 section (query) onto the Day0 section (reference) by
    affine + B-spline free-form deformation guided by the embeddings, then matching in the
    REGISTERED frame. The earlier `cast` wrapper skipped CAST_STACK and chose candidates in raw
    slide coordinates (different slides -> arbitrary), so it must not be reported as CAST.
    Coupling is SOFT: each Day9 cell spreads unit mass over its k nearest Day0 cells in the
    registered frame, weighted by a Gaussian kernel of registered distance (scale = median
    nearest-neighbour distance) times the non-negative embedding cosine similarity. No integer
    ties. Smoke-test on one unit before submitting the array: reg_params fields and the
    gpu/device convention are taken from the CAST source and may need adjusting."""
    import torch
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import normalize
    from CAST import CAST_MARK, CAST_STACK
    try:
        from CAST import reg_params
    except ImportError:
        from CAST.CAST_Stack import reg_params
    k = int(p.get("knn_phys", 20))
    A, B, _, _ = bc.load_pair(p)
    A_xy = np.asarray(A.obsm["spatial"], dtype=np.float64); B_xy = np.asarray(B.obsm["spatial"], dtype=np.float64)
    # CAST_MARK takes plain arrays: coords -> delaunay_dgl(np.array(df)), expression -> torch.tensor(..., device).
    # Demo1 feeds normalize_total(1e4) counts; adata.X here is log-normalised, which CAST also accepts.
    coords = {"A": A_xy, "B": B_xy}
    expr = {"A": bc._todense(A.X).astype(np.float32), "B": bc._todense(B.X).astype(np.float32)}
    tag = str(p.get("_urepr") or p.get("slice_id", "FIRST"))
    tag = "".join(c if (c.isalnum() or c in "._-") else "_" for c in tag)
    work = os.path.join(p["out_dir"], "cast_stack", tag); os.makedirs(work, exist_ok=True)
    embed = CAST_MARK(coords, expr, work, task_name_t=tag, if_plot=False)      # dict {"A": tensor, "B": tensor}
    # --- registration: query B (Day9) onto reference A (Day0). Parameters exactly as CAST's demo2
    # (demo2_CAST_Stack_Align_S4_to_S1): a bare reg_params(dataname=...) is NOT valid because the
    # list fields (d_list, attention_params, iterations_bs, ...) default to [] and CAST_STACK indexes them.
    params = reg_params(dataname=tag,
                        gpu=0 if torch.cuda.is_available() else -1,
                        diff_step=5,
                        iterations=500, dist_penalty1=0, bleeding=500,
                        d_list=[3, 2, 1, 1 / 2, 1 / 3], attention_params=[None, 3, 1, 0],
                        dist_penalty2=[0], alpha_basis_bs=[500], meshsize=[8], iterations_bs=[400],
                        attention_params_bs=[[None, 3, 1, 0]], mesh_weight=[None])
    params.alpha_basis = torch.Tensor([1 / 1000, 1 / 1000, 1 / 50, 5, 5]).reshape(5, 1).to(params.device)
    coords_raw = {"A": A_xy.astype(np.float64), "B": B_xy.astype(np.float64)}            # numpy, as CAST_STACK expects (.mean(0), .shape)
    embed_cpu = {kk: v.detach().cpu() for kk, v in embed.items()}                          # CAST_STACK calls .cpu() on these
    coords_final = CAST_STACK(coords_raw, embed_cpu, work, graph_list=["B", "A"], params_dist=params, mid_visual=False)
    B_reg = np.asarray(coords_final["B"].numpy() if hasattr(coords_final["B"], "numpy") else coords_final["B"], dtype=np.float64)   # query, mapped into A's frame
    A_ref = np.asarray(coords_final["A"], dtype=np.float64)                                                                        # reference, unchanged (rescale factor 1)
    # --- matching in the registered frame ---
    embA = normalize(embed["A"].detach().cpu().numpy()); embB = normalize(embed["B"].detach().cpu().numpy())
    nn = NearestNeighbors(n_neighbors=min(k, A.n_obs)).fit(A_ref)
    dist, cand = nn.kneighbors(B_reg)
    scale = np.median(dist[:, 0]) + 1e-9
    P = np.zeros((A.n_obs, B.n_obs), dtype=np.float32)
    for j in range(B.n_obs):
        w = np.exp(-0.5 * (dist[j] / scale) ** 2) * np.clip(embB[j] @ embA[cand[j]].T, 0, None)
        if w.sum() <= 0: w = np.exp(-0.5 * (dist[j] / scale) ** 2)
        P[cand[j], j] += (w / w.sum()).astype(np.float32)
    return bc.make_typematrix(P, A.obs["Tier2"].astype(str).values,
                              B.obs["Tier2"].astype(str).values,
                              display_types=bc.type_set(p),
                              row_normalize=p.get("row_normalize", False))


METHODS = {"moscot": m_moscot, "paste2": m_paste2,
           "stalign": m_stalign, "cast": m_cast,
           "moscot_ub": m_moscot_ub, "slat": m_slat, "cast_stack": m_cast_stack}


def _run_one_unit(method_name, p):
    """Worker entry point for the parallel unit loop. Returns
    (unit_repr, matrix_or_None, error_str_or_None) so failures survive as a
    pooled `None` instead of killing the whole run (same semantics as serial)."""
    urepr = p.get("_urepr", "?")
    try:
        m = METHODS[method_name](p)
        gc.collect()
        return (urepr, m, None)
    except Exception as e:
        gc.collect()
        return (urepr, None, f"{type(e).__name__}: {e}")


def _partial_path(out_dir, method_name, lineage, i):
    """Where one unit's RAW per-pair matrix is written (design-agnostic)."""
    return os.path.join(out_dir, "partials",
                        f"{method_name}_{lineage}_unit{i:03d}.csv")


def _save_partial(out_dir, method_name, lineage, i, M):
    """Persist one unit's raw per-pair matrix (output #1: the untouched output)."""
    if M is None or getattr(M, "size", 0) == 0:
        return
    os.makedirs(os.path.join(out_dir, "partials"), exist_ok=True)
    M.to_csv(_partial_path(out_dir, method_name, lineage, i))


# --------------------------------------------------------------------------- #
# Concordance Index (output #3). Ground truth = ordinal lineage geometry; the
# metric is the fraction of biology-defined orderings the similarity respects.
# It reads ONLY within-row order, so it is identical whether a matrix is raw or
# row-normalized -> all methods (raw here) score on one footing.
# --------------------------------------------------------------------------- #
def _gt_distance(a, b, axis, L):
    """Ground-truth dissimilarity: 0 if same type; |pos-diff| for two on-axis
    types; L (> any on-axis gap) if either is off-axis."""
    if a == b:
        return 0
    if a in axis and b in axis:
        return abs(axis[a] - axis[b])
    return L


def _concordance_matrix(M, axis, L):
    """One slice's concordance tally + per-comparison ledger rows.
    Ledger row = (query, closer_type, farther_type, is_identity, outcome)
    where outcome in {1 concordant, 0 discordant, 0.5 tie}."""
    conc = comp = id_conc = id_comp = gr_conc = gr_comp = 0.0
    ledger = []
    for a in list(M.index):
        cands = [c for c in M.columns if np.isfinite(M.loc[a, c])]
        for x in range(len(cands)):
            for y in range(x + 1, len(cands)):
                b, c = cands[x], cands[y]
                gb, gc = _gt_distance(a, b, axis, L), _gt_distance(a, c, axis, L)
                if gb == gc:
                    continue                       # tie in ground truth -> not comparable
                closer, farther = (b, c) if gb < gc else (c, b)
                sb, sc = M.loc[a, closer], M.loc[a, farther]
                out = 1.0 if sb > sc else (0.0 if sb < sc else 0.5)
                is_id = (min(gb, gc) == 0)         # closer is the self type
                comp += 1; conc += out
                if is_id:
                    id_comp += 1; id_conc += out
                else:
                    gr_comp += 1; gr_conc += out
                ledger.append((a, closer, farther, is_id, out))
    return dict(conc=conc, comp=comp, id_conc=id_conc, id_comp=id_comp,
                gr_conc=gr_conc, gr_comp=gr_comp, ledger=ledger)


def _ci(conc, comp):
    return (conc / comp) if comp else float("nan")


def _write_ci_outputs(mats, axis, method_name, lineage, out_dir):
    """Write the CI bundle (output #3) from a list of per-slice raw matrices:
      <method>_ci_summary_<lineage>.csv  -> per-slice CI + identity/gradient + overall
      <method>_ci_ledger_<lineage>.csv   -> each biology-defined ordering, correct in N/M slices
      <method>_rowdiag_<lineage>.csv     -> per query type: self-sim, self-rank, top-1 %
    """
    mats = [m for m in mats if m is not None and getattr(m, "size", 0)]
    if not axis or not mats:
        print(f"  [ci/{lineage}] no axis or no matrices; skip", flush=True)
        return
    L = len(axis)                                  # > any on-axis distance (= len-1)

    # --- per-slice CI + identity/gradient split ---
    rows, ledger_acc = [], {}
    for si, M in enumerate(mats):
        t = _concordance_matrix(M, axis, L)
        rows.append(dict(slice_index=si,
                         CI=_ci(t["conc"], t["comp"]),
                         CI_identity=_ci(t["id_conc"], t["id_comp"]),
                         CI_gradient=_ci(t["gr_conc"], t["gr_comp"]),
                         n_comparisons=int(t["comp"])))
        for (q, closer, farther, is_id, out) in t["ledger"]:
            k = (q, closer, farther, is_id)
            acc = ledger_acc.setdefault(k, [0.0, 0])
            acc[0] += out; acc[1] += 1
    summ = pd.DataFrame(rows)
    overall = dict(slice_index="ALL(mean)",
                   CI=summ["CI"].mean(),
                   CI_identity=summ["CI_identity"].mean(),
                   CI_gradient=summ["CI_gradient"].mean(),
                   n_comparisons=int(summ["n_comparisons"].sum()))
    summ = pd.concat([summ, pd.DataFrame([overall])], ignore_index=True)
    summ.to_csv(os.path.join(out_dir, f"{method_name}_ci_summary_{lineage}.csv"),
                index=False)

    # --- comparison ledger: which orderings hold, across slices ---
    led = [dict(query=q, more_similar_to=closer, than_to=farther,
                relation=("identity" if is_id else "gradient"),
                correct=round(s, 2), n_slices=n,
                frac=round(s / n, 3) if n else float("nan"))
           for (q, closer, farther, is_id), (s, n) in ledger_acc.items()]
    (pd.DataFrame(led).sort_values(["relation", "frac", "query"])
       .to_csv(os.path.join(out_dir, f"{method_name}_ci_ledger_{lineage}.csv"),
               index=False))

    # --- per-type row diagnostics ---
    qtypes = sorted(set().union(*[set(m.index) for m in mats]))
    diag = []
    for q in qtypes:
        selfsims, selfranks, top1 = [], [], []
        argmaxes = {}
        for M in mats:
            if q not in M.index:
                continue
            row = M.loc[q].dropna()
            if not len(row):
                continue
            if q in row.index:
                selfsims.append(float(row[q]))
                # rank of self (1 = largest); ties -> min rank
                selfranks.append(int((row > row[q]).sum()) + 1)
                top1.append(1.0 if row.idxmax() == q else 0.0)
            argmaxes[row.idxmax()] = argmaxes.get(row.idxmax(), 0) + 1
        if not selfsims:
            continue
        diag.append(dict(query=q,
                         on_axis=(q in axis),
                         mean_self_sim=round(float(np.mean(selfsims)), 4),
                         mean_self_rank=round(float(np.mean(selfranks)), 2),
                         top1_pct=round(100 * float(np.mean(top1))),
                         most_common_match=max(argmaxes, key=argmaxes.get),
                         n_slices=len(selfsims)))
    pd.DataFrame(diag).to_csv(
        os.path.join(out_dir, f"{method_name}_rowdiag_{lineage}.csv"), index=False)
    ov = summ.iloc[-1]
    print(f"  [ci/{lineage}] CI={ov['CI']:.3f} "
          f"(identity={ov['CI_identity']:.3f}, gradient={ov['CI_gradient']:.3f})",
          flush=True)


def run_method(method_name, base, only_unit=None, pool_only=False):
    """Run one method over both lineages for base['design'], pooling all units and
    saving <method>_typematrix_<lineage>.csv into base['out_dir']. Returns
    {lineage: DataFrame}. Callable from an Rmd (via reticulate) or from main().

    base keys: design, adata_path, sample_type, seed, max_cells_per_half, out_dir,
    lineages (list), and method knobs (raster_dx, niter, sigma_mult, knn_phys).

    Array-job modes (for slow methods like PASTE2, one SLURM task per pair):
      only_unit=N  -> run ONLY unit N (both lineages), write its per-unit matrix
                      to partials/<method>_<lineage>_unitNNN.csv, do NOT pool.
      pool_only=True -> read all partials for each lineage and write the final
                      pooled <method>_typematrix_<lineage>.csv.
    Default (both None/False) is unchanged: in-process loop + pool + save.
    """
    method = METHODS[method_name]
    design = str(base.get("design", "halfsplit"))
    os.makedirs(base["out_dir"], exist_ok=True)
    results = {}
    for lineage in list(base.get("lineages", ["Fibroblast", "Epithelial"])):
        p0 = dict(base); p0["lineage"] = lineage; p0.setdefault("slice_id", "FIRST")
        units = bc.all_units(p0)   # slice ids (halfsplit) or 45 slice pairs (slicepair)

        # ---- array-job modes: run one unit, or pool the per-unit partials ----
        if only_unit is not None or pool_only:
            os.makedirs(os.path.join(base["out_dir"], "partials"), exist_ok=True)
            if pool_only:
                mats = []
                for i in range(len(units)):
                    f = _partial_path(base["out_dir"], method_name, lineage, i)
                    if os.path.exists(f):
                        mats.append(pd.read_csv(f, index_col=0))
                    else:
                        print(f"  [pool/{lineage}] MISSING unit {i} ({units[i]})",
                              flush=True)
                        mats.append(None)
                have = sum(m is not None for m in mats)
                print(f"  [pool/{lineage}] pooling {have}/{len(units)} partials",
                      flush=True)
                # output #2: simple pooled average (no normalization)
                pooled = bc.pool_matrices(mats, display_types=bc.type_set(p0))
                bc.save_typematrix(pooled, base["out_dir"],
                                   f"{method_name}_typematrix_{lineage}.csv",
                                   row_normalize=False)
                # output #3: Concordance Index bundle from the per-slice partials
                # (skip for day0day9 — cross-condition uses the downstream rank-recovery
                #  scorer, not the within-condition ordinal CI)
                if design != "day0day9":
                    _write_ci_outputs(mats, bc.lineage_axis(p0), method_name,
                                      lineage, base["out_dir"])
                results[lineage] = pooled
            else:
                i = int(only_unit)
                if i < 0 or i >= len(units):
                    print(f"  [unit {i}/{lineage}] out of range 0..{len(units)-1}; "
                          f"skip", flush=True)
                    continue
                u = units[i]
                p = dict(p0)
                if design in ("slicepair", "day0day9"):
                    p["slice_a"], p["slice_b"] = u
                else:
                    p["slice_id"] = u
                p["_urepr"] = str(u)   # unique tag (CAST uses it for its mark dir)
                f = _partial_path(base["out_dir"], method_name, lineage, i)
                try:
                    M = method(p)
                    M.to_csv(f)
                    print(f"  [unit {i}/{lineage}] {u}: ok -> {f}", flush=True)
                except Exception as e:
                    print(f"  [unit {i}/{lineage}] {u} FAILED: "
                          f"{type(e).__name__}: {e}", flush=True)
                gc.collect()
            continue   # array modes never fall through to the in-process loop

        print(f"\n=== [{method_name}/{design}] {lineage}: pooling {len(units)} units ===",
              flush=True)
        # Build the per-unit parameter dicts up front.
        unit_ps = []
        for u in units:
            p = dict(p0)
            if design in ("slicepair", "day0day9"):
                p["slice_a"], p["slice_b"] = u
            else:
                p["slice_id"] = u
            p["_urepr"] = str(u)
            unit_ps.append(p)

        # Pooling averages over the list, so unit order does not matter -> the
        # loop can run in parallel. Opt in with BENCH_NJOBS>1 (default 1 keeps
        # the exact serial behaviour for every other method/design).
        njobs = int(os.environ.get("BENCH_NJOBS", "1"))
        mats = []
        if njobs > 1 and len(unit_ps) > 1:
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor, as_completed
            print(f"  [parallel] {len(unit_ps)} units across {njobs} workers",
                  flush=True)
            ctx = mp.get_context("spawn")  # light top-level imports -> cheap workers
            with ProcessPoolExecutor(max_workers=njobs, mp_context=ctx) as ex:
                # map each future back to its unit index so partials are named right
                futs = {ex.submit(_run_one_unit, method_name, p): idx
                        for idx, p in enumerate(unit_ps)}
                for fut in as_completed(futs):
                    idx = futs[fut]
                    urepr, mat, err = fut.result()
                    mats.append(mat)
                    _save_partial(base["out_dir"], method_name, lineage, idx, mat)
                    print(f"  {urepr}: {'ok' if err is None else 'FAILED: ' + err}",
                          flush=True)
        else:
            for idx, p in enumerate(unit_ps):
                try:
                    M = method(p)
                    mats.append(M)
                    _save_partial(base["out_dir"], method_name, lineage, idx, M)
                    print(f"  {p['_urepr']}: ok", flush=True)
                except Exception as e:
                    print(f"  {p['_urepr']} FAILED: {type(e).__name__}: {e}",
                          flush=True)
                    mats.append(None)
                gc.collect()
        # output #1 (raw per-pair matrices) is written above into partials/.
        # output #2: simple pooled average (no normalization).
        pooled = bc.pool_matrices(mats, display_types=bc.type_set(p0))
        bc.save_typematrix(pooled, base["out_dir"],
                           f"{method_name}_typematrix_{lineage}.csv", row_normalize=False)
        # output #3: Concordance Index bundle (skip for cross-condition day0day9).
        if design != "day0day9":
            _write_ci_outputs(mats, bc.lineage_axis(p0), method_name, lineage,
                              base["out_dir"])
        results[lineage] = pooled
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--method", required=True, choices=list(METHODS))
    ap.add_argument("--adata", required=True, help="path to adata.h5ad")
    ap.add_argument("--out", required=True, help="output dir for the CSVs")
    ap.add_argument("--max_cells", type=int, default=4000, help="cells per half (cap)")
    ap.add_argument("--sample_type", default="Healthy", help="Day0 = Healthy")
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--lineages", nargs="+", default=["Fibroblast", "Epithelial"])
    # STalign / CAST knobs (ignored by the others)
    ap.add_argument("--raster_dx", type=float, default=30)
    ap.add_argument("--niter", type=int, default=200)
    ap.add_argument("--sigma_mult", type=float, default=2)
    ap.add_argument("--knn_phys", type=int, default=50)
    ap.add_argument("--design", default="halfsplit",
                    choices=["halfsplit", "slicepair", "day0day9"],
                    help="halfsplit = Benchmark A (one slice -> halves); "
                         "slicepair = Benchmark B (all C(n,2) full-slice pairs); "
                         "day0day9 = cross-condition Day0 x Day9 grid (Fibroblast)")
    ap.add_argument("--sample_type_a", default="Healthy", help="day0day9: A condition (Day0)")
    ap.add_argument("--sample_type_b", default="DSS9", help="day0day9: B condition (Day9)")
    ap.add_argument("--day0day9_min_cells", type=int, default=30,
                    help="day0day9: keep Tier2 types with >= this many cells per side per pair")
    ap.add_argument("--http_proxy", default=os.environ.get("http_proxy", ""))
    ap.add_argument("--unit", type=int, default=None,
                    help="array mode: run ONLY this unit index (e.g. SLURM_ARRAY_TASK_ID) "
                         "and write its partial matrix; do not pool")
    ap.add_argument("--pool", action="store_true",
                    help="array mode: pool the per-unit partials into the final CSV")
    ap.add_argument("--row_normalize", action="store_true",
                    help="row-normalize each per-pair matrix (default OFF = raw output, "
                         "so the per-pair CSVs and the pooled average are untouched)")
    args = ap.parse_args()

    bc.configure_proxy(dict(http_proxy=args.http_proxy, https_proxy=args.http_proxy,
                            no_proxy="localhost,127.0.0.1"))
    run_method(args.method, dict(
        adata_path=args.adata, sample_type=args.sample_type, design=args.design,
        seed=args.seed, max_cells_per_half=args.max_cells, out_dir=args.out,
        lineages=list(args.lineages), raster_dx=args.raster_dx, niter=args.niter,
        sigma_mult=args.sigma_mult, knn_phys=args.knn_phys,
        row_normalize=args.row_normalize,
        sample_type_a=args.sample_type_a, sample_type_b=args.sample_type_b,
        day0day9_min_cells=args.day0day9_min_cells),
        only_unit=args.unit, pool_only=args.pool)


if __name__ == "__main__":
    main()
