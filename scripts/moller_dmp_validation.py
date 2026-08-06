#!/usr/bin/env python3
"""Moller 2025 562-cluster-specific DMP probe-level validation across 3 independent EPIC cohorts.

Loads the EXACT 562 DMPs (cg IDs + cluster membership k1-k5 + methylation trend) extracted
from Møller 2025 Supplementary_File_S2.xlsx, and validates at the probe level whether those
DMPs recapitulate Møller's methylation trends and discriminate the high-risk (k4-like) group
in our three cohorts: GSE283928 (tumor clusters), GSE207937 (NS vs secreting), Neou 2020 (internal clusters).
"""
import os, json, gzip, csv, time, math
import numpy as np
from collections import Counter

BASE = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
PROJ = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
OUT = os.path.join(PROJ, "results", "Moller2025")
os.makedirs(OUT, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. Load the 562 DMPs (Moller 2025 Supp File S2, "Specific probe methylation")
# ---------------------------------------------------------------------------
import openpyxl, warnings
warnings.filterwarnings("ignore")
S2 = os.path.join("G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Moller2025_supplement", "extracted", "Supplementary_File_S2.xlsx")
wb = openpyxl.load_workbook(S2, read_only=True)
ws = wb["Specific probe methylation"]
dmp = []
for i, r in enumerate(ws.iter_rows(values_only=True)):
    if i == 0:
        continue
    dmp.append({"probe": r[1], "cluster": r[2], "family": r[3],
                "enhancer": r[4], "cpg_loc": r[5], "trend": r[6]})
print(f"[1] Loaded {len(dmp)} Moller DMPs. Memberships={dict(Counter(d['cluster'] for d in dmp))} "
      f"Trend={dict(Counter(d['trend'] for d in dmp))}")
DMP_PROBES = [d["probe"] for d in dmp]
# predicted sign from trend: Hyper -> group mean higher (+1); Hypo -> lower (-1)
SIGN = np.array([1.0 if d["trend"] == "Hypermethylated" else -1.0 for d in dmp])
CLUS = np.array([d["cluster"] for d in dmp])
# k4 (highest-risk, HR 3.63) specific DMPs
K4MASK = (CLUS == "k4")

# ---------------------------------------------------------------------------
# 2. Loaders (reuse exact logic from prior scripts)
# ---------------------------------------------------------------------------
def load_neou():
    bdir = os.path.join(BASE, "E-MTAB-7762", "beta_out")
    binpath = os.path.join(bdir, "neou_beta_full.bin")
    with open(binpath, "rb") as f:
        nrow = int(np.frombuffer(f.read(4), dtype="<i4")[0])
        ncol = int(np.frombuffer(f.read(4), dtype="<i4")[0])
        arr = np.frombuffer(f.read(), dtype="<f8").reshape(nrow, ncol)   # probes x samples (R col-major)
    X = arr.T.astype(np.float32)   # samples x probes
    probes = [l.strip() for l in open(os.path.join(bdir, "neou_beta_probes.txt"))]
    samples = [l.strip() for l in open(os.path.join(bdir, "neou_beta_samples.txt"))]
    colmean = np.nanmean(X, axis=0, keepdims=True)
    X = np.where(np.isnan(X), colmean, X).astype(np.float32)
    print(f"    Neou: {X.shape[0]} samples x {X.shape[1]} probes")
    return np.array(probes), np.array(samples), X  # probes, samples, (n_samp x n_probe)

def load_gse283928():
    MAT = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
    rows, probes, header = [], [], None
    with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!") or line.startswith("#") or not line.strip():
                continue
            p = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if header is None:
                header = p; continue
            probes.append(p[0]); rows.append(p[1:])
    X = np.array(rows, dtype=np.float32)              # probes x cols
    hdr = header[1:]
    beta_cols = [j for j, h in enumerate(hdr) if h.endswith("_Beta")]
    dpv_cols = [j for j, h in enumerate(hdr) if "Detection Pval" in h]
    print(f"    GSE283928: cols={len(hdr)} beta={len(beta_cols)} dpval={len(dpv_cols)}")
    Xb = X[:, beta_cols]                              # probes x 118 (beta only)
    beta_ids = [hdr[j] for j in beta_cols]             # keep full e.g. "Case1_Beta"
    # labels: only _Beta rows from cluster_assignment (exclude detection-pval artifact -> cluster 1)
    clu = {}
    with open(os.path.join(BASE, "results", "GSE283928", "cluster_assignment_k5.csv")) as f:
        for row in csv.DictReader(f):
            clu[row["sample_id"]] = int(row["cluster"])
    lab = np.array([clu.get(bid, -1) for bid in beta_ids])   # -1 for missing
    keep = lab >= 0
    Xb = Xb[:, keep]; beta_ids = [beta_ids[j] for j in range(len(beta_ids)) if keep[j]]
    lab = lab[keep]
    print(f"    tumor samples kept={Xb.shape[1]} clusters={dict(Counter(lab.tolist()))}")
    return np.array(probes), np.array(beta_ids), Xb, lab   # probes, samples, (probe x n), labels

def load_gse207937():
    MAT = os.path.join(BASE, "GSE207937_series_matrix.txt.gz")
    rows, probes, header = [], [], None
    with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!") or line.startswith("#") or not line.strip():
                continue
            p = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if header is None:
                header = p; continue
            probes.append(p[0]); rows.append(p[1:])
    X = np.array(rows, dtype=np.float32)
    Xb = (1.0 / (1.0 + np.power(2.0, -X))).astype(np.float32)   # M -> beta
    del X
    sample_ids = np.array([h.replace('"', '') for h in header[1:]])
    # nfpt labels from sample_labels.csv
    nfpt = {}
    with open(os.path.join(PROJ, "results", "GSE207937", "sample_labels.csv")) as f:
        for row in csv.DictReader(f):
            nfpt[row["sample_id"]] = int(row["nfpt"])
    lab = np.array([nfpt.get(s, -1) for s in sample_ids])
    keep = lab >= 0
    Xb = Xb[:, keep]; sample_ids = sample_ids[keep]; lab = lab[keep]
    print(f"    GSE207937: samples={Xb.shape[1]} nfpt=1:{int((lab==1).sum())} nfpt=0:{int((lab==0).sum())}")
    return np.array(probes), sample_ids, Xb, lab

# ---------------------------------------------------------------------------
# 3. Build per-cohort DMP beta submatrix (562 x n_samp)
# ---------------------------------------------------------------------------
def dmp_submatrix(probes_arr, X_probe_x_samp):
    # X_probe_x_samp: (n_probe x n_samp). Return full (562 x n_samp) matrix, NaN where absent.
    pidx = {p: i for i, p in enumerate(probes_arr)}
    out = np.full((len(DMP_PROBES), X_probe_x_samp.shape[1]), np.nan)
    present = np.zeros(len(DMP_PROBES), dtype=bool)
    for j, p in enumerate(DMP_PROBES):
        i = pidx.get(p, -1)
        if i >= 0:
            out[j] = X_probe_x_samp[i, :]
            present[j] = True
    return out, present

cohorts = {}
# Neou
pN, sN, XN = load_neou()
subN, presN = dmp_submatrix(pN, XN.T)   # XN is samp x probe -> transpose
# Neou high-risk = internal_cluster 0
ic = {}
with open(os.path.join(PROJ, "results", "Neou", "neou_internal_clusters.csv")) as f:
    for row in csv.DictReader(f):
        ic[row["sample_id"]] = int(row["internal_cluster"])
H_N = np.array([ic.get(s, -1) == 0 for s in sN])
cohorts["Neou"] = {"beta": subN, "present": presN, "H": H_N, "n": subN.shape[1]}

# GSE283928
p2, s2, X2, lab2 = load_gse283928()
keep2 = np.isin(lab2, [0, 3, 4])     # drop c2 (n=2 artifact)
X2k = X2[:, keep2]; lab2k = lab2[keep2]
sub2, pres2 = dmp_submatrix(p2, X2k)
H_2 = (lab2k == 0)   # c0 = high-risk (k4 analog) vs SF1 clusters c3/c4
cohorts["GSE283928"] = {"beta": sub2, "present": pres2, "H": H_2, "n": sub2.shape[1]}

# GSE207937
p7, s7, X7, lab7 = load_gse207937()
sub7, pres7 = dmp_submatrix(p7, X7)
H_7 = (lab7 == 1)   # NS = NFPT = high-risk analog
cohorts["GSE207937"] = {"beta": sub7, "present": pres7, "H": H_7, "n": sub7.shape[1]}

# union present across cohorts
present_union = np.zeros(len(DMP_PROBES), dtype=bool)
for c in cohorts.values():
    present_union |= c["present"]
print(f"[3] DMPs present in all 3 cohorts: {int(np.all(np.array([c['present'] for c in cohorts.values()]),axis=0).sum())}/562; "
      f"present in >=1: {int(present_union.sum())}/562")

# ---------------------------------------------------------------------------
# 4. Per-cohort directionality + effect (Welch t, normal approx)
# ---------------------------------------------------------------------------
def norm_cdf(x): return 0.5 * (1 + math.erf(x / math.sqrt(2)))

def welch_z(a, b):
    ma, mb = a.mean(), b.mean()
    va, vb = a.var(ddof=1), b.var(ddof=1)
    na, nb = len(a), len(b)
    se = math.sqrt(va/na + vb/nb) + 1e-12
    t = (ma - mb) / se
    return t, ma - mb

def roc_auc(score, y):
    pos = score[y]; neg = score[~y]
    n1, n2 = len(pos), len(neg)
    if n1 == 0 or n2 == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    order = np.argsort(allv, kind="mergesort")
    ranks = np.empty(len(allv)); ranks[order] = np.arange(1, len(allv)+1)
    R1 = ranks[:n1].sum()
    return (R1 - n1*(n1+1)/2) / (n1*n2)

results = {"cohorts": {}, "per_probe": {}}
meta_z = np.zeros((len(DMP_PROBES), 3))   # signed z per probe per cohort
meta_valid = np.zeros((len(DMP_PROBES), 3), dtype=bool)
cohort_names = ["Neou", "GSE283928", "GSE207937"]

for ci, cname in enumerate(cohort_names):
    C = cohorts[cname]
    beta = C["beta"]; H = C["H"]; pres = C["present"]
    nh = int(H.sum()); nl = int((~H).sum())
    delta = np.full(len(DMP_PROBES), np.nan)
    z = np.full(len(DMP_PROBES), np.nan)
    correct = np.full(len(DMP_PROBES), np.nan)
    for j in range(len(DMP_PROBES)):
        if not pres[j]:
            continue
        v = beta[j]
        a = v[H]; b = v[~H]
        if len(a) < 3 or len(b) < 3:
            continue
        t, d = welch_z(a, b)
        delta[j] = d
        s = SIGN[j]
        # signed z: positive (one-sided p/2) if delta in predicted direction, negative if opposite
        sz = abs(norm_cdf(-abs(t))) if d * s >= 0 else -abs(norm_cdf(-abs(t)))
        meta_z[j, ci] = sz
        meta_valid[j, ci] = True
        correct[j] = 1.0 if d * s > 0 else 0.0
    # direction consistency (binomial vs 0.5) over present probes
    cvals = correct[~np.isnan(correct)]
    ncorr = int(cvals.sum()); ntot = len(cvals)
    p_binom = sum(math.comb(ntot, k) * 0.5**ntot for k in range(ncorr, ntot+1)) if ntot > 0 else 1.0
    # k4 subset
    cvals_k4 = correct[K4MASK & pres]
    ncorr_k4 = int(cvals_k4.sum()); ntot_k4 = len(cvals_k4)
    p_binom_k4 = sum(math.comb(ntot_k4, k) * 0.5**ntot_k4 for k in range(ncorr_k4, ntot_k4+1)) if ntot_k4 > 0 else 1.0
    # DMP signature score AUC (hyper - hypo)
    hyper = SIGN > 0
    score = np.full(beta.shape[1], np.nan)
    for s in range(beta.shape[1]):
        score[s] = np.nanmean(beta[hyper & pres, s]) - np.nanmean(beta[~hyper & pres, s])
    auc = roc_auc(score, H)
    results["cohorts"][cname] = {
        "n_samples": C["n"], "n_high_risk": nh, "n_low_risk": nl,
        "n_dmp_present": ntot, "frac_correct_direction": ncorr/ntot if ntot else None,
        "binom_p_all": p_binom, "n_correct_all": ncorr,
        "n_dmp_present_k4": ntot_k4, "frac_correct_k4": ncorr_k4/ntot_k4 if ntot_k4 else None,
        "binom_p_k4": p_binom_k4, "n_correct_k4": ncorr_k4,
        "signature_auc": auc,
    }
    print(f"    {cname}: present={ntot} frac_correct={ncorr/ntot:.3f} (binom p={p_binom:.2e}) "
          f"| k4: {ncorr_k4}/{ntot_k4} (p={p_binom_k4:.2e}) | signature AUC={auc:.3f}")

# ---------------------------------------------------------------------------
# 5. Cross-cohort Stouffer meta per probe + BH FDR
# ---------------------------------------------------------------------------
n_present_probe = meta_valid.sum(axis=1) > 0
z_meta = np.full(len(DMP_PROBES), np.nan)
p_meta = np.full(len(DMP_PROBES), np.nan)
for j in range(len(DMP_PROBES)):
    vals = meta_z[j][meta_valid[j]]
    if len(vals) == 0:
        continue
    k = len(vals)
    zm = vals.sum() / math.sqrt(k)
    z_meta[j] = zm
    p_meta[j] = 2 * (1 - norm_cdf(abs(zm)))
# BH FDR
mask = ~np.isnan(p_meta)
pvals = p_meta[mask]
order = np.argsort(pvals)
m = len(pvals)
fdr = np.empty(m)
prev = 1.0
for i in range(m-1, -1, -1):
    r = pvals[order[i]] * m / (i+1)
    prev = min(prev, r)
    fdr[i] = prev
fdr_sorted = np.empty(m); fdr_sorted[order] = fdr
fdr_all = np.full(len(DMP_PROBES), np.nan); fdr_all[mask] = fdr_sorted
sig = (p_meta < 0.05) & (meta_valid.sum(axis=1) >= 2)
n_sig = int(sig.sum())
n_sig_correct = int((sig & (SIGN * np.sign(z_meta) > 0)).sum())
print(f"[5] Cross-cohort meta: probes with meta-p<0.05 (>=2 cohorts) = {n_sig}/{int(mask.sum())}; "
      f"of those, direction-consistent with Moller = {n_sig_correct}")

# ---------------------------------------------------------------------------
# 6. Save outputs
# ---------------------------------------------------------------------------
# per-probe table
with open(os.path.join(OUT, "moller_dmp_probe_validation.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["probe", "moller_cluster", "trend", "predicted_sign",
                "present_Neou", "present_GSE283928", "present_GSE207937",
                "delta_Neou", "delta_GSE283928", "delta_GSE207937",
                "z_meta", "p_meta", "fdr", "direction_consistent"])
    for j in range(len(DMP_PROBES)):
        w.writerow([DMP_PROBES[j], CLUS[j], dmp[j]["trend"], int(SIGN[j]),
                    int(cohorts["Neou"]["present"][j]), int(cohorts["GSE283928"]["present"][j]),
                    int(cohorts["GSE207937"]["present"][j]),
                    "" if np.isnan(meta_z[j,0]) else round(meta_z[j,0],4),
                    "" if np.isnan(meta_z[j,1]) else round(meta_z[j,1],4),
                    "" if np.isnan(meta_z[j,2]) else round(meta_z[j,2],4),
                    "" if np.isnan(z_meta[j]) else round(z_meta[j],4),
                    "" if np.isnan(p_meta[j]) else f"{p_meta[j]:.3e}",
                    "" if np.isnan(fdr_all[j]) else round(fdr_all[j],4),
                    int(SIGN[j]*np.sign(z_meta[j])>0) if not np.isnan(z_meta[j]) else ""])

summary = {
    "n_dmp": len(DMP_PROBES),
    "memberships": dict(Counter(CLUS.tolist())),
    "trend": dict(Counter(dmp[j]["trend"] for j in range(len(dmp)))),
    "present_union": int(present_union.sum()),
    "present_all3": int(np.all(np.array([c["present"] for c in cohorts.values()]), axis=0).sum()),
    "cohorts": results["cohorts"],
    "meta_n_sig_p05": n_sig,
    "meta_n_sig_direction_correct": n_sig_correct,
}
json.dump(summary, open(os.path.join(OUT, "moller_dmp_validation_summary.json"), "w"), indent=2)
print("[6] Saved moller_dmp_probe_validation.csv + summary.json")

# stash for figure
np.save(os.path.join(OUT, "_meta_z.npy"), meta_z)
np.save(os.path.join(OUT, "_z_meta.npy"), z_meta)
np.save(os.path.join(OUT, "_p_meta.npy"), p_meta)
np.save(os.path.join(OUT, "_fdr.npy"), fdr_all)
json.dump({"SIGN": SIGN.tolist(), "CLUS": CLUS.tolist(), "DMP_PROBES": DMP_PROBES},
          open(os.path.join(OUT, "_dmp_meta.json"), "w"))
print("DONE")
