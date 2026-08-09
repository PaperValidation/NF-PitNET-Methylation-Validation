#!/usr/bin/env python3
"""Definitive C<->k mapping using Moller 2025's 562 cluster-specific DMPs across all 3 cohorts.
For each Moller cluster k (k1-k5) and each cohort group, compute fraction of k's DMPs whose
methylation direction matches Moller's labeled trend (binomial test). Best-matching group = k's analog.
This resolves whether GSE283928 c0 is the high-risk k4 or an SF1-lineage cluster.
"""
import os, csv, gzip, json, math
import numpy as np
from collections import Counter

from _paths import DATASET_PATH as BASE, DATASET_PATH_TEMP as TEMP
PROJ = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
OUT = os.path.join(PROJ, "results", "Moller2025")

import openpyxl, warnings
warnings.filterwarnings("ignore")
wb = openpyxl.load_workbook(os.path.join("G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Moller2025_supplement", "extracted", "Supplementary_File_S2.xlsx"), read_only=True)
ws = wb["Specific probe methylation"]
dmp = []
for i, r in enumerate(ws.iter_rows(values_only=True)):
    if i == 0: continue
    dmp.append({"probe": r[1], "cluster": r[2], "trend": r[6]})
DMP_PROBES = [d["probe"] for d in dmp]
SIGN = np.array([1.0 if d["trend"] == "Hypermethylated" else -1.0 for d in dmp])
CLUS = np.array([d["cluster"] for d in dmp])
MOLLER_CLUSTERS = ["k1", "k2", "k3", "k4", "k5"]

# ---- loaders (reuse) ----
def load_neou():
    bdir = os.path.join(BASE, "E-MTAB-7762", "beta_out")
    with open(os.path.join(bdir, "neou_beta_full.bin"), "rb") as f:
        nrow = int(np.frombuffer(f.read(4), dtype="<i4")[0])
        ncol = int(np.frombuffer(f.read(4), dtype="<i4")[0])
        arr = np.frombuffer(f.read(), dtype="<f8").reshape(nrow, ncol)
    X = arr.T.astype(np.float32)
    probes = [l.strip() for l in open(os.path.join(bdir, "neou_beta_probes.txt"))]
    samples = [l.strip() for l in open(os.path.join(bdir, "neou_beta_samples.txt"))]
    cm = np.nanmean(X, axis=0, keepdims=True); X = np.where(np.isnan(X), cm, X)
    return np.array(probes), np.array(samples), X

def load_gse283928():
    MAT = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
    rows, probes, header = [], [], None
    with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!") or line.startswith("#") or not line.strip(): continue
            p = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if header is None: header = p; continue
            probes.append(p[0]); rows.append(p[1:])
    X = np.array(rows, dtype=np.float32); hdr = header[1:]
    bc = [j for j, h in enumerate(hdr) if h.endswith("_Beta")]
    Xb = X[:, bc]; bids = [hdr[j] for j in bc]
    clu = {}
    with open(os.path.join(TEMP, "results", "GSE283928", "cluster_assignment_k5.csv")) as f:
        for row in csv.DictReader(f): clu[row["sample_id"]] = int(row["cluster"])
    lab = np.array([clu.get(b, -1) for b in bids]); keep = np.isin(lab, [0,2,3,4])
    return np.array(probes), Xb[:, keep], lab[keep]

def load_gse207937():
    MAT = os.path.join(BASE, "GSE207937_series_matrix.txt.gz")
    rows, probes, header = [], [], None
    with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!") or line.startswith("#") or not line.strip(): continue
            p = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if header is None: header = p; continue
            probes.append(p[0]); rows.append(p[1:])
    X = np.array(rows, dtype=np.float32)
    Xb = (1.0/(1.0+np.power(2.0,-X))).astype(np.float32); del X
    sids = np.array([h.replace('"','') for h in header[1:]])
    nfpt = {}
    with open(os.path.join(PROJ, "results", "GSE207937", "sample_labels.csv")) as f:
        for row in csv.DictReader(f): nfpt[row["sample_id"]] = int(row["nfpt"])
    lab = np.array([nfpt.get(s, -1) for s in sids]); keep = lab >= 0
    return np.array(probes), Xb[:, keep], lab[keep]

def dmp_matrix(probes_arr, X):
    pidx = {p: i for i, p in enumerate(probes_arr)}
    out = np.full((len(DMP_PROBES), X.shape[1]), np.nan)
    for j, p in enumerate(DMP_PROBES):
        i = pidx.get(p, -1)
        if i >= 0: out[j] = X[i]
    return out

pN, sN, XN = load_neou(); DN = dmp_matrix(pN, XN.T)
p2, X2, lab2 = load_gse283928(); D2 = dmp_matrix(p2, X2)
p7, X7, lab7 = load_gse207937(); D7 = dmp_matrix(p7, X7)

def binom_p(correct, tot):
    return sum(math.comb(tot, m) * 0.5**tot for m in range(correct, tot+1)) if tot else 1.0

def group_match(D, labels, groups):
    """For each Moller cluster k, find best-matching group + frac_correct + p."""
    res = {}
    for k in MOLLER_CLUSTERS:
        kmask = (CLUS == k) & ~np.isnan(D[:, 0])
        nk = int(kmask.sum())
        best = None
        per = {}
        for g in groups:
            H = (labels == g); L = ~H
            if H.sum() < 3 or L.sum() < 3: continue
            corr = tot = 0
            for j in np.where(kmask)[0]:
                a = D[j][H]; b = D[j][L]
                if len(a) < 3 or len(b) < 3: continue
                if (a.mean()-b.mean()) * SIGN[j] > 0: corr += 1
                tot += 1
            frac = corr/tot if tot else 0
            p = binom_p(corr, tot)
            per[str(g)] = {"frac": round(frac,3), "n": tot, "p": f"{p:.1e}", "correct": corr}
            # raw best by frac (record regardless of significance)
            if best is None or frac > best[1]:
                best = (g, frac, p, corr, tot)
        sig = (best is not None and best[1] > 0.5 and best[2] < 0.05)
        res[k] = {"nDMP": nk, "per_group": per,
                  "best_match": (best[0], round(best[1],3), best[2], best[3], best[4]) if best else None,
                  "significant": bool(sig)}
    return res

# GSE283928 groups: c0,c2,c3,c4 ; Neou: int0,int1,int2 ; GSE207937: ns(1), secreting(0)
m283 = group_match(D2, lab2, [0, 2, 3, 4])

ic = {}
with open(os.path.join(PROJ, "results", "Neou", "neou_internal_clusters.csv")) as f:
    for row in csv.DictReader(f):
        ic[row["sample_id"]] = int(row["internal_cluster"])
labN = np.array([ic.get(s, -1) for s in sN])
mNeu = group_match(DN, labN, [0, 1, 2])

lab7g = np.where(lab7 == 1, "NS", "secreting")
m207 = group_match(D7, lab7g, ["NS", "secreting"])

# ---- report ----
print("\n================ MOLLER k -> COHORT GROUP MAPPING ================")
for name, mm in [("GSE283928(c0/c2/c3/c4)", m283), ("Neou(int0/int1/int2)", mNeu),
                 ("GSE207937(NS/secreting)", m207)]:
    print(f"\n### {name}")
    for k in MOLLER_CLUSTERS:
        r = mm[k]
        bm = r["best_match"]
        if bm is None:
            bm_s = "NONE"
        else:
            bm_s = f"c{bm[0]}" if isinstance(bm[0], int) else str(bm[0])
        sig = " *SIG*" if r.get("significant") else ""
        print(f"  {k} (nDMP={r['nDMP']}): best match = {bm_s} "
              f"(frac={bm[1] if bm else '-'}, p={bm[2] if bm else '-'}){sig} | per-group {r['per_group']}")

# ---- k4 signature transferability (high-risk k4 = TPIT/PIT1) ----
print("\n================ k4 (HR 3.63, highest-risk) DMP SIGNATURE TRANSFERABILITY ================")
k4mask = (CLUS == "k4")
for name, D, labels, grps, glab in [
    ("GSE283928", D2, lab2, [0,2,3,4], lambda g: f"c{g}"),
    ("Neou", DN, labN, [0,1,2], lambda g: f"int{g}"),
    ("GSE207937", D7, lab7g, ["NS","secreting"], lambda g: str(g))]:
    hyper = (SIGN > 0) & k4mask; hypo = (SIGN < 0) & k4mask
    print(f"  {name}: mean(k4-hyper)-mean(k4-hypo) signature per group:")
    for g in grps:
        if name == "GSE207937":
            m = (labels == g)
        else:
            m = (labels == g)
        if m.sum() < 3: continue
        sc = np.nanmean(D[hyper][:, m]) - np.nanmean(D[hypo][:, m])
        print(f"     {glab(g)}: k4-signature={sc:+.4f}  n={int(m.sum())}")

# ---- save ----
summary = {"GSE283928": m283, "Neou": mNeu, "GSE207937": m207}
json.dump(summary, open(os.path.join(OUT, "moller_k_cluster_mapping.json"), "w"), indent=2, default=str)
print("\nSaved moller_k_cluster_mapping.json")

