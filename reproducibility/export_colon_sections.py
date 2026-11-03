#!/usr/bin/env python3
"""Export the per-section cell tables used by every colon stage.
    python3 export_colon_sections.py [adata.h5ad] [all_slices]
Reads the obs table of the atlas h5ad (Cadinu et al. 2024, Dryad doi:10.5061/dryad.rjdfn2zh3) and writes one CSV per
Slice_ID (52 sections) with the columns the pipeline expects: x, y, Mouse_ID, Technical_repeat_number, Sample_type,
Slice_ID, FOV, Cell_ID, sample, N_genes, Tier1, Tier2, Tier3, Leiden_neigh and the UMAP columns of the atlas. Rows keep
the atlas order; the first (unnamed) column is the atlas obs index. Also writes atlas_sections.csv (one row per
section: Slice_ID, Mouse_ID, Technical_repeat_number, Sample_type, n_cells) next to the output folder. Only h5py and pandas are needed."""
import sys, os, h5py, numpy as np, pandas as pd
H5 = sys.argv[1] if len(sys.argv) > 1 else "adata.h5ad"; OUT = sys.argv[2] if len(sys.argv) > 2 else "all_slices"
COLS = ["x", "y", "Mouse_ID", "Technical_repeat_number", "Sample_type", "Slice_ID", "FOV", "Cell_ID", "sample", "N_genes",
        "Tier1", "Tier2", "Tier3", "Leiden_neigh", "Neigh_umap_x", "Neigh_umap_y", "Tier2_umap_x", "Tier2_umap_y",
        "Tier3_umap_x", "Tier3_umap_y", "Tier1_umap_x", "Tier1_umap_y"]
def column(obs, k):
    g = obs[k]
    if isinstance(g, h5py.Group):                      # categorical
        cats = g["categories"][:]; cats = cats.astype(str) if cats.dtype.kind in "SO" else cats
        return pd.Categorical.from_codes(g["codes"][:], cats).astype(object if cats.dtype.kind in "OU" else cats.dtype)
    v = g[:]; return v.astype(str) if v.dtype.kind in "SO" else v
with h5py.File(H5) as h:
    obs = h["obs"]; idx = obs["_index"][:].astype(str)
    df = pd.DataFrame({k: column(obs, k) for k in COLS}, index=idx)
os.makedirs(OUT, exist_ok=True)
for sid, d in df.groupby("Slice_ID", sort=False, observed=True):
    d.to_csv(os.path.join(OUT, f"{sid}.csv"))
sec = (df.groupby("Slice_ID", sort=True, observed=True).agg(Mouse_ID=("Mouse_ID", "first"), Technical_repeat_number=("Technical_repeat_number", "first"),
       Sample_type=("Sample_type", "first"), n_cells=("Cell_ID", "size")).reset_index())
sec.to_csv(os.path.join(os.path.dirname(os.path.abspath(OUT)), "atlas_sections.csv"), index=False)   # section table used by the figure scripts
print(f"wrote {df['Slice_ID'].nunique()} section tables ({len(df):,} cells) to {OUT}/ and atlas_sections.csv")
