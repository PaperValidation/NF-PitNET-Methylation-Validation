#!/usr/bin/env python3
# ============================================================================
# Build & cache the GSE283928 (Belakhoua 2025) cluster centroids used for the
# Møller-2025 strict projection test.
# ----------------------------------------------------------------------------
# Reads the big GSE283928 normalized EPIC matrix ONCE, computes per-cluster
# centroids (c0/c2/c3/c4; cluster 1 = detection-pvalue artifact, EXCLUDED) and
# per-probe mean/SD over the 118 real tumor samples, and caches everything to a
# small .npz so the projection step never reloads the 917 MB matrix.
#
# Output: G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Moller2025_derived/gse283928_centroids.npz
#   - probes      : (P,) cg probe ids (intersection-ready with Møller cg)
#   - centroid_c0/c2/c3/c4 : (P,) mean beta per probe for each real cluster
#   - probe_mean  : (P,) probe mean across real tumor samples
#   - probe_std   : (P,) probe SD   across real tumor samples
#
# Run: python scripts/build_gse283928_centroid.py
# ============================================================================
import os, time, numpy as np, pandas as pd

from _paths import DATASET_PATH as BASE_G
MAT    = os.path.join(BASE_G, "GSE283928_MatrixProcessed_Normalized.txt.gz")
CL     = os.path.join(BASE_G, "results", "GSE283928", "cluster_assignment_k5.csv")
OUT    = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Moller2025_derived/gse283928_centroids.npz"

# cluster id -> name; 1 is the detection-pvalue artifact column group -> drop
REAL = {0: "c0", 2: "c2", 3: "c3", 4: "c4"}

print("[1] Loading GSE283928 matrix (this is ~917 MB gz, may take 1-2 min)...")
t0 = __import__("time").time()
# NOTE: do NOT pass a scalar dtype here — it would also coerce the ProbeName
# index column to float and crash. Read, then cast data to float32.
mat = pd.read_csv(MAT, sep="\t", index_col=0, compression="gzip")
mat = mat.astype(np.float32)
print(f"    matrix {mat.shape}, {time.time()-t0:.1f}s")

cl = pd.read_csv(CL)  # columns: sample_id, cluster
cl = cl[cl["cluster"].isin(REAL)]          # drop artifact cluster 1
cl["name"] = cl["cluster"].map(REAL)
present = set(cl["sample_id"])
mat_cols = [c for c in mat.columns if c in present]
mat = mat[mat_cols]
cl = cl.set_index("sample_id").loc[mat_cols].reset_index()
print(f"    real tumor samples kept: {mat.shape[1]}")

# keep only cg probes (match Møller's cg set)
cg_mask = mat.index.str.startswith("cg")
probes = np.array([str(s) for s in mat.index[cg_mask]], dtype="U20")  # clean unicode, no pickle
M = mat.values[cg_mask]                    # (P, N) float32
print(f"    cg probes: {M.shape[0]}, samples: {M.shape[1]}")

# per-probe scale from real tumor samples
probe_mean = M.mean(axis=1)
probe_std  = M.std(axis=1)
probe_std[probe_std == 0] = 1.0

centroids = {}
for cid, name in REAL.items():
    cols = cl.loc[cl["cluster"] == cid, "sample_id"].tolist()
    centroids[f"centroid_{name}"] = M[:, [mat.columns.get_loc(c) for c in cols]].mean(axis=1)

np.savez_compressed(
    OUT,
    probes=probes,
    probe_mean=probe_mean.astype(np.float32),
    probe_std=probe_std.astype(np.float32),
    **{k: v.astype(np.float32) for k, v in centroids.items()},
)
print(f"[2] Saved centroids -> {OUT}")
print("    centroid sample sizes:", {REAL[c]: int((cl['cluster']==c).sum()) for c in REAL})
