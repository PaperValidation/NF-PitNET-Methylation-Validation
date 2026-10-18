#!/usr/bin/env python3
"""Diagnostic: map each Moller cluster k1-k5 to the best-matching GSE283928 tumor cluster
using the cluster-specific DMP directions. Resolves whether c0 == k4."""
import os, csv, gzip, json
import numpy as np
from collections import Counter

from _paths import DATASET_PATH as BASE, DATASET_PATH_TEMP as TEMP
PROJ = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
OUT = os.path.join(PROJ, "results", "Moller2025")

import openpyxl, warnings
warnings.filterwarnings("ignore")
wb = openpyxl.load_workbook(os.path.join(BASE, "Moller2025_supplement", "extracted", "Supplementary_File_S2.xlsx"), read_only=True)
ws = wb["Specific probe methylation"]
dmp = []
for i, r in enumerate(ws.iter_rows(values_only=True)):
    if i == 0: continue
    dmp.append({"probe": r[1], "cluster": r[2], "trend": r[6]})
DMP_PROBES = [d["probe"] for d in dmp]
SIGN = np.array([1.0 if d["trend"] == "Hypermethylated" else -1.0 for d in dmp])
CLUS = np.array([d["cluster"] for d in dmp])

# load GSE283928 beta at DMPs + labels
MAT = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
rows, probes, header = [], [], None
with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!") or line.startswith("#") or not line.strip(): continue
        p = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
        if header is None: header = p; continue
        probes.append(p[0]); rows.append(p[1:])
X = np.array(rows, dtype=np.float32)
hdr = header[1:]
beta_cols = [j for j, h in enumerate(hdr) if h.endswith("_Beta")]
Xb = X[:, beta_cols]
beta_ids = [hdr[j] for j in beta_cols]
clu = {}
with open(os.path.join(TEMP, "results", "GSE283928", "cluster_assignment_k5.csv")) as f:
    for row in csv.DictReader(f):
        clu[row["sample_id"]] = int(row["cluster"])
lab = np.array([clu.get(b, -1) for b in beta_ids])
keep = np.isin(lab, [0, 2, 3, 4])
Xb = Xb[:, keep]; lab = lab[keep]; beta_ids = [beta_ids[j] for j in range(len(beta_ids)) if keep[j]]
pidx = {p: i for i, p in enumerate(probes)}
D = np.full((len(DMP_PROBES), Xb.shape[1]), np.nan)
for j, p in enumerate(DMP_PROBES):
    if p in pidx: D[j] = Xb[pidx[p]]

gse_clusters = [0, 2, 3, 4]
print("Moller cluster k -> best-matching GSE283928 cluster (by frac_correct of k's DMPs):")
mapping = {}
for k in ["k1", "k2", "k3", "k4", "k5"]:
    kmask = (CLUS == k) & ~np.isnan(D[:, 0])
    nk = int(kmask.sum())
    best_c, best_frac, best_p = None, -1, 1
    detail = {}
    for c in gse_clusters:
        H = (lab == c); L = ~H
        if H.sum() < 3 or L.sum() < 3: continue
        correct = 0; tot = 0
        for j in np.where(kmask)[0]:
            a = D[j][H]; b = D[j][L]
            if len(a) < 3 or len(b) < 3: continue
            d = a.mean() - b.mean()
            if d * SIGN[j] > 0: correct += 1
            tot += 1
        frac = correct / tot if tot else 0
        # binomial p
        pval = sum(__import__("math").comb(tot, m) * 0.5**tot for m in range(correct, tot+1)) if tot else 1
        detail[c] = (round(frac, 3), tot, f"{pval:.1e}")
        if frac > best_frac:
            best_frac, best_c, best_p = frac, c, pval
    mapping[k] = (best_c, round(best_frac, 3))
    print(f"  {k}: nDMP={nk}  per-GSE283928-cluster frac_correct={detail}  -> best match c{best_c} (frac={best_frac:.3f})")

# Also: signature score per GSE283928 cluster for k4 DMPs specifically
print("\nGSE283928 cluster mean(k4-hyper beta) - mean(k4-hypo beta) [k4 signature; higher=more k4-like]:")
k4mask = (CLUS == "k4")
hyper = (SIGN > 0) & k4mask
hypo = (SIGN < 0) & k4mask
for c in gse_clusters:
    m = lab == c
    sc = np.nanmean(D[hyper][:, m]) - np.nanmean(D[hypo][:, m]) if (hyper.any() and hypo.any()) else float("nan")
    # also mean raw beta at k4 DMPs
    raw = np.nanmean(D[k4mask][:, m])
    print(f"  c{c}: k4-signature={sc:.4f}  mean_raw_beta_at_k4DMPs={raw:.4f}  n={int(m.sum())}")

# Møller-key gene beta per GSE283928 cluster (recompute for cross-check)
print("\nCross-check: mean beta at Møller-key genes per GSE283928 cluster (from prior: c0 highest):")
# reuse lineage_marker_probes.csv Moller_key set
lin = os.path.join(TEMP, "results", "GSE283928", "lineage_marker_probes.csv")
mk = set()
with open(lin) as f:
    for row in csv.DictReader(f):
        if row.get("set") == "Moller_key" or row.get("Set") == "Moller_key":
            mk.add(row["probe"])
mkidx = [pidx[p] for p in mk if p in pidx]
for c in gse_clusters:
    m = lab == c
    print(f"  c{c}: mean Moller_key beta={np.mean(Xb[mkidx][:, m]):.4f}  n={int(m.sum())}")
