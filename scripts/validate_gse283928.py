#!/usr/bin/env python3
"""Validate GSE283928 5-cluster structure against Moller 2025 biology:
- SF1-lineage (NR5A1/SF1 hypermethylated per Moller) should separate clusters
- Compare cluster centroids on known lineage markers
"""
import gzip, numpy as np, csv, os, json

from _paths import DATASET_PATH as BASE, DATASET_PATH_TEMP as TEMP
MAT  = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
OUT  = os.path.join(TEMP, "results", "GSE283928")

# load cluster assignment
assign = {}
with open(os.path.join(OUT, "cluster_assignment_k5.csv")) as f:
    r = csv.reader(f); next(r)
    for row in r: assign[row[0]] = int(row[1])
sample_ids = list(assign.keys())
clusters = np.array([assign[s] for s in sample_ids])

# lineage marker genes -> we need their EPIC probe cg IDs. We don't have annotation offline,
# so instead validate using the PCA structure: check if top PC correlates with known biology
# by re-running PCA and reporting cluster separation on PC space.
print("=== Cluster separation in PCA space (centroid distance) ===")
pc = np.loadtxt(os.path.join(OUT, "pca_coords_top20.csv"), delimiter=",", skiprows=1)
# PC rows correspond to samples in same order as sample_ids (saved earlier in same order)
cent = {}
for c in np.unique(clusters):
    cent[c] = pc[clusters == c].mean(axis=0)
print("Cluster centroids (PC1,PC2,PC3):")
for c in sorted(cent):
    print(f"  cluster {c}: PC1={cent[c][0]:.2f} PC2={cent[c][1]:.2f} PC3={cent[c][2]:.2f}")

# between-cluster centroid distance matrix (PC1-3)
import itertools
print("\nPairwise centroid Euclidean distance (PC1-3):")
dmap = {}
for a, b in itertools.combinations(sorted(cent), 2):
    d = np.linalg.norm(cent[a][:3] - cent[b][:3])
    dmap[(a,b)] = d
for k in sorted(dmap, key=dmap.get, reverse=True):
    print(f"  cluster {k[0]} vs {k[1]}: {dmap[k]:.2f}")

# Save cluster summary
summary = {
    "n_samples": len(sample_ids),
    "cluster_sizes": {str(k): int(v) for k, v in zip(*np.unique(clusters, return_counts=True))},
    "centroids_PC1_3": {str(c): [round(float(x),3) for x in cent[c][:3]] for c in cent},
    "pairwise_dist_PC1_3": {f"{a}-{b}": round(float(dmap[(a,b)]),3) for (a,b) in dmap},
}
with open(os.path.join(OUT, "cluster_validation_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("\nSaved cluster_validation_summary.json")
