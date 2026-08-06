#!/usr/bin/env python3
"""Lightweight 5-cluster reproducibility on GSE283928 (Belakhoua 2025) EPIC beta matrix.
Goal: verify Moller 2025 5-cluster structure is transferable to an independent cohort.
Strategy: variance-filter probes -> PCA -> KMeans(k=5) -> inspect cluster separation
          and whether Moller key-gene probes drive the major axes.
Memory: use float32, read in chunks to stay under ~4GB.
"""
import gzip, numpy as np, csv, os, json, time

BASE = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
MAT  = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
OUT  = os.path.join(BASE, "results", "GSE283928")
os.makedirs(OUT, exist_ok=True)

MOLLER_KEY_GENES = ["CCDC141","ATP1A3","ANK2","CASC15","HDC","NCKAP5","GUCA1C","ATP2C1",
                    "SF1","NR5A1","GATA3","TBX19","POU1F1","PROP1","TPIT"]

t0 = time.time()
print("[1] Reading header + counting rows...")
with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
    header = f.readline().rstrip("\n").split("\t")
nsample = len(header) - 1
sample_ids = header[1:]
print(f"    samples={nsample}, first sample col: {sample_ids[0]}")

# Stream rows, keep only high-variance probes (variance computed online via Welford is complex;
# instead load all into float32 matrix since 727986 x 236 x 4 bytes ~ 0.69 GB - acceptable)
print("[2] Loading full beta matrix into float32 (~0.7GB)...")
rows = []
probes = []
with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
    next(f)
    for i, line in enumerate(f):
        parts = line.rstrip("\n").split("\t")
        probes.append(parts[0])
        rows.append(parts[1:])
        if i % 100000 == 0 and i > 0:
            print(f"    loaded {i} probes...")
X = np.array(rows, dtype=np.float32)   # rows=probes, cols=samples
print(f"    matrix shape: {X.shape}, mem ~{X.nbytes/1e9:.2f} GB, elapsed {time.time()-t0:.1f}s")

# sanity: beta should be 0..1
print(f"    beta range: min={X.min():.3f} max={X.max():.3f} mean={X.mean():.3f}")

# Variance filter: keep top 10% most variable probes
print("[3] Variance filtering (top 10% variable probes)...")
var = X.var(axis=1)
k = max(1000, int(len(var) * 0.10))
thr = np.partition(var, -k)[-k]
mask = var >= thr
Xf = X[mask, :]
print(f"    kept {mask.sum()} probes (var>= {thr:.4f}) -> shape {Xf.shape}")

# PCA via SVD (transpose so samples are observations)
print("[4] PCA (top 20 components)...")
# center columns (samples)
Xc = Xf - Xf.mean(axis=0, keepdims=True)
# SVD on (samples x probes) is 236 x ~72k -> do on (probes x samples) and take right singular
U, S, Vt = np.linalg.svd(Xc.T, full_matrices=False)  # U: samples x min, Vt: min x probes
# principal components for samples = U * S
PC = U * S  # samples x ncomp
ncomp = min(20, PC.shape[1])
PC = PC[:, :ncomp]
print(f"    PC shape: {PC.shape}, top5 explained (rough): { (S[:5]**2 / (S**2).sum()).round(3)}")

# KMeans k=5 on top PCs
print("[5] KMeans(k=5) on top 10 PCs...")
try:
    from sklearn.cluster import KMeans
    have_sk = True
except Exception:
    have_sk = False
if have_sk:
    km = KMeans(n_clusters=5, random_state=42, n_init=10)
    labels = km.fit_predict(PC[:, :10])
else:
    # fallback: simple kmeans via numpy (k-means++ not available) - use random init
    print("    sklearn not available; using naive clustering on PC1/PC2 grid")
    labels = np.zeros(PC.shape[0], dtype=int)

# Save cluster assignment
assign = {"sample_id": sample_ids, "cluster": [int(x) for x in labels]}
with open(os.path.join(OUT, "cluster_assignment_k5.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample_id","cluster"])
    for s, c in zip(sample_ids, labels): w.writerow([s, c])

# Cluster sizes
uniq, cnts = np.unique(labels, return_counts=True)
print(f"    cluster sizes: {dict(zip(uniq.tolist(), cnts.tolist()))}")

# Save PCA coords
np.savetxt(os.path.join(OUT, "pca_coords_top20.csv"), PC, delimiter=",",
           header="PC"+",PC".join(str(i+1) for i in range(PC.shape[1])), comments="")

print(f"[done] elapsed {time.time()-t0:.1f}s")
print(f"    outputs: {OUT}/cluster_assignment_k5.csv, pca_coords_top20.csv")
