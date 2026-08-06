#!/usr/bin/env python3
"""Neou 2020 (E-MTAB-7762) cross-cohort methylation validation.

Validates the MOLLER 2025 / BELAKHOUA 2025 EPIC methylation taxonomy on Neou 2020,
a THIRD independent EPIC 850K cohort (correcting the prior SDRF mislabel that called it 450K).
Methods mirror gse207937_analysis.py for direct comparability:
  - joint 3-cohort PCA (Neou + GSE283928 tumor + GSE207937) in common-probe beta space
  - nearest-centroid projection of Neou (and GSE207937) onto GSE283928 4 tumor clusters
  - lineage (SF1/TPIT/PIT1) + Moller-2025 key-gene methylation scoring on Neou NFPT vs functioning
NOTE: GSE283928 local matrix = 118 tumor Beta cols + 118 Detection-Pval cols (artifact).
      The 118 tumor cols are used; the 4 tumor clusters are {0:70,2:2,3:22,4:24}.
"""
import gzip, csv, json, time, re, os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

BASE = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
NE = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
MAT_283 = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
MAT_207 = os.path.join(BASE, "GSE207937_series_matrix.txt.gz")
SAMP_207 = os.path.join(BASE, "samples", "GSE207937_samples.csv")
CLUST_283 = os.path.join(BASE, "results", "GSE283928", "cluster_assignment_k5.csv")
LIN = os.path.join(BASE, "results", "GSE283928", "lineage_marker_probes.csv")
NEU_ANNO = os.path.join(BASE, "E-MTAB-7762", "beta_out", "neou_sample_annotation.csv")
NEU_BIN = os.path.join(BASE, "E-MTAB-7762", "beta_out", "neou_beta_full.bin")
NEU_PROBES = os.path.join(BASE, "E-MTAB-7762", "beta_out", "neou_beta_probes.txt")
NEU_SAMPLES = os.path.join(BASE, "E-MTAB-7762", "beta_out", "neou_beta_samples.txt")
OUT = os.path.join(NE, "results", "Neou")
os.makedirs(OUT, exist_ok=True)

def clean(s):
    return s.strip().strip('"')

def m2beta(m):
    return 1.0 / (1.0 + np.power(2.0, -m))

def impute(M):
    """Replace NaN with column (probe) mean; fallback 0.5 if an entire probe is NaN."""
    col_mean = np.nanmean(M, axis=0)
    col_mean = np.where(np.isnan(col_mean), 0.5, col_mean)
    idx = np.where(np.isnan(M))
    if len(idx[0]):
        M = M.copy()
        M[idx] = np.take(col_mean, idx[1])
    return M

def perm_p(a, b, n=20000, seed=0):
    """two-sided permutation p-value for difference in means (scipy-free)."""
    rng = np.random.default_rng(seed)
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    obs = a.mean() - b.mean()
    pool = np.concatenate([a, b])
    cnt = 0
    for _ in range(n):
        rng.shuffle(pool)
        d = pool[:len(a)].mean() - pool[len(a):].mean()
        if abs(d) >= abs(obs):
            cnt += 1
    return (cnt + 1) / (n + 1)

t0 = time.time()

# ---------------------------------------------------------------------------
# 1. Neou beta (samples x probes) from R binary export
# ---------------------------------------------------------------------------
print("[1] Loading Neou beta binary ...")
with open(NEU_BIN, "rb") as f:
    nrow = int(np.fromfile(f, dtype=np.int32, count=1)[0])
    ncol = int(np.fromfile(f, dtype=np.int32, count=1)[0])
    data = np.fromfile(f, dtype=np.float64, count=nrow * ncol)
Bmat = data.reshape(nrow, ncol, order="F")      # probes x samples (R column-major)
Neou = np.ascontiguousarray(Bmat.T)             # samples x probes
with open(NEU_PROBES, encoding="utf-8") as f:
    neu_probes = [l.rstrip("\n") for l in f]
with open(NEU_SAMPLES, encoding="utf-8") as f:
    neu_samples = [l.rstrip("\n") for l in f]
neu_p2i = {p: i for i, p in enumerate(neu_probes)}
print(f"    Neou shape {Neou.shape}; probes={len(neu_probes)} samples={len(neu_samples)}")
Neou = impute(Neou)
print(f"    elapsed {time.time()-t0:.1f}s")

# Neou annotation
print("[2] Neou sample annotation ...")
neu_anno = {}
with open(NEU_ANNO, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        neu_anno[row["individual"]] = row
neu_nfpt = np.array([int(neu_anno.get(s, {}).get("nfpt", 0)) for s in neu_samples])
neu_stage = [neu_anno.get(s, {}).get("disease_staging", "NA") for s in neu_samples]
neu_lineage = [neu_anno.get(s, {}).get("lineage", "NA") for s in neu_samples]
print(f"    NFPT={int(neu_nfpt.sum())}/{len(neu_samples)} ; stages:",
      {k: neu_stage.count(k) for k in set(neu_stage)})

# ---------------------------------------------------------------------------
# 3. GSE283928 (Beta cols only)
# ---------------------------------------------------------------------------
print("[3] Loading GSE283928 (Beta cols) ...")
with gzip.open(MAT_283, "rt", encoding="utf-8", errors="ignore") as f:
    hdr = f.readline().rstrip("\n").split("\t")
beta_idx = [i for i, h in enumerate(hdr) if h.endswith("_Beta")]
beta_cols = [hdr[i] for i in beta_idx]
print(f"    {len(beta_cols)} tumor Beta cols")
rows = []; probes283 = []
with gzip.open(MAT_283, "rt", encoding="utf-8", errors="ignore") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        probes283.append(p[0])
        rows.append([p[i] for i in beta_idx])
X283 = np.array(rows, dtype=np.float32).T        # samples(118) x probes
X283 = impute(X283)
p283 = np.array(probes283)
del rows
p283_i = {p: j for j, p in enumerate(p283)}
print(f"    GSE283928 shape {X283.shape}; elapsed {time.time()-t0:.1f}s")

# GSE283928 tumor cluster labels (filter to _Beta rows only)
print("[4] GSE283928 tumor cluster labels ...")
clu = {}
with open(CLUST_283, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["sample_id"].endswith("_Beta"):
            clu[row["sample_id"]] = int(row["cluster"])
lab283 = np.array([clu.get(c, -1) for c in beta_cols])
tumor_clusters = sorted(set(lab283.tolist()) - {-1})
print(f"    tumor clusters {tumor_clusters}; sizes",
      {int(c): int((lab283 == c).sum()) for c in tumor_clusters})

# ---------------------------------------------------------------------------
# 5. GSE207937 (M -> beta) + subtype labels
# ---------------------------------------------------------------------------
print("[5] Loading GSE207937 (M->beta) ...")
h2 = None; p207 = []; r207 = []
with gzip.open(MAT_207, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!") or line.startswith("#") or not line.strip():
            continue
        parts = [clean(x) for x in line.rstrip("\n").split("\t")]
        if h2 is None:
            h2 = parts; continue
        p207.append(parts[0]); r207.append(parts[1:])
X207m = np.array(r207, dtype=np.float32)          # probes x 77 (M-values)
X207 = m2beta(X207m).astype(np.float32)           # beta
X207 = X207.T                                      # samples(77) x probes
X207 = impute(X207)
p207_arr = np.array(p207)
p207_i = {p: j for j, p in enumerate(p207)}
gse207_samples = h2[1:]
print(f"    GSE207937 shape {X207.shape}")
sub_by_gsm = {}
with open(SAMP_207, encoding="utf-8", errors="ignore") as f:
    for row in csv.DictReader(f):
        gsm = clean(row.get("sample_geo_accession", ""))
        ch = row.get("sample_characteristics_ch1", "") or ""
        m = re.search(r"subtypes:\s*([A-Za-z\-]+)", ch)
        sub_by_gsm[gsm] = m.group(1) if m else "NA"
gse207_sub = [sub_by_gsm.get(s, "NA") for s in gse207_samples]
gse207_nfpt = np.array([1 if s == "NS" else 0 for s in gse207_sub])
print(f"    GSE207937 subtypes:", {k: gse207_sub.count(k) for k in set(gse207_sub)})

# ---------------------------------------------------------------------------
# 6. Common probes across all three; build combined matrix
# ---------------------------------------------------------------------------
print("[6] Intersect probes + combined matrix ...")
set_n = set(neu_probes); set_2 = set(p283.tolist()); set_7 = set(p207)
common = sorted(set_n & set_2 & set_7)
ci_n = [neu_p2i[p] for p in common]
ci_2 = [p283_i[p] for p in common]
ci_7 = [p207_i[p] for p in common]
Nc = Neou[:, ci_n]          # 86 x P
X2c = X283[:, ci_2]         # 118 x P
X7c = X207[:, ci_7]         # 77 x P
combined = np.vstack([Nc, X2c, X7c]).astype(np.float32)   # 281 x P
print(f"    common probes = {len(common)}; combined {combined.shape}; elapsed {time.time()-t0:.1f}s")

# ---------------------------------------------------------------------------
# 7. Variance filter + joint PCA
# ---------------------------------------------------------------------------
print("[7] Variance filter (top 10%) + PCA(20) ...")
var = combined.var(axis=0)
k = max(1000, int(len(var) * 0.10))
thr = np.partition(var, -k)[-k]
mask = var >= thr
Xf = combined[:, mask]
print(f"    kept {int(mask.sum())} probes")
pca = PCA(n_components=20, random_state=42)
PC = pca.fit_transform(Xf)      # 281 x 20
evr = pca.explained_variance_ratio_
print(f"    PC1/PC2 explained {evr[0]:.3f}/{evr[1]:.3f}; cum20={evr[:20].sum():.3f}")

# cohort tags for combined rows
cohort = (["Neou"] * Nc.shape[0]) + (["GSE283928"] * X2c.shape[0]) + (["GSE207937"] * X7c.shape[0])
label_nfpt = np.concatenate([neu_nfpt, np.full(X2c.shape[0], -1), gse207_nfpt])  # -1 = no NFPT label (GSE283928)
sample_ids = neu_samples + beta_cols + gse207_samples
# save joint PCA coords
with open(os.path.join(OUT, "joint_pca_coords.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["cohort", "sample_id", "nfpt"] + [f"PC{i+1}" for i in range(20)])
    for i in range(PC.shape[0]):
        w.writerow([cohort[i], sample_ids[i], int(label_nfpt[i])] + [f"{x:.5f}" for x in PC[i]])

# ---------------------------------------------------------------------------
# 8. Cross-cohort nearest-centroid projection onto GSE283928 tumor clusters
# ---------------------------------------------------------------------------
print("[8] Cross-cohort projection onto GSE283928 tumor centroids ...")
# centroids in common-probe beta space
Cmat = np.zeros((len(tumor_clusters), len(common)), dtype=np.float32)
for j, c in enumerate(tumor_clusters):
    m = lab283 == c
    Cmat[j] = X2c[m, :].mean(axis=0)
cn = np.linalg.norm(Cmat, axis=1)
def assign(Xq):
    D = np.zeros((Xq.shape[0], len(tumor_clusters)), dtype=np.float32)
    for i in range(Xq.shape[0]):
        D[i] = np.linalg.norm(Cmat - Xq[i:i+1], axis=1)
    a_e = np.array(tumor_clusters)[np.argmin(D, axis=1)]
    xn = np.linalg.norm(Xq, axis=1)
    corr = (Xq @ Cmat.T) / (xn[:, None] * cn[None, :])
    a_c = np.array(tumor_clusters)[np.argmax(corr, axis=1)]
    return a_e, a_c, D
neu_e, neu_c, neu_D = assign(Nc)
g207_e, g207_c, g207_D = assign(X7c)

# Neou assignment table
with open(os.path.join(OUT, "neou_crosscohort_assignment.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample_id", "disease_staging", "nfpt", "lineage",
                "assigned_cluster_eucl", "assigned_cluster_corr"] +
               [f"dist_c{c}" for c in tumor_clusters])
    for i, s in enumerate(neu_samples):
        a = neu_anno.get(s, {})
        w.writerow([s, a.get("disease_staging", ""), neu_nfpt[i], a.get("lineage", ""),
                    int(neu_e[i]), int(neu_c[i])] + [f"{x:.4f}" for x in neu_D[i]])

# GSE207937 per-sample assignment (for heatmap figure)
with open(os.path.join(OUT, "gse207937_crosscohort_assignment.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample_id", "subtype", "nfpt", "assigned_cluster_eucl",
                "assigned_cluster_corr"] + [f"dist_c{c}" for c in tumor_clusters])
    for i, s in enumerate(gse207_samples):
        w.writerow([s, gse207_sub[i], int(gse207_nfpt[i]), int(g207_e[i]),
                    int(g207_c[i])] + [f"{x:.4f}" for x in g207_D[i]])

# contingency: Neou nfpt vs assigned cluster
def contingency(assigned, nfpt_flag):
    out = {}
    for st, fl in [("NFPT", 1), ("Functioning", 0)]:
        m = nfpt_flag == fl
        vals, cc = np.unique(assigned[m], return_counts=True)
        out[st] = {int(v): int(c) for v, c in zip(vals, cc)}
    return out
neu_ct = contingency(neu_e, neu_nfpt)
# 3-way: cohort x assigned cluster
threeway = {}
for name, ass, n in [("Neou", neu_e, neu_nfpt), ("GSE207937", g207_e, gse207_nfpt)]:
    threeway[name] = {int(c): int((ass == c).sum()) for c in tumor_clusters}
print("    Neou NFPT vs cluster:", neu_ct)
print("    3-way cluster counts:", threeway)

# ---------------------------------------------------------------------------
# 9. Lineage + Moller-key methylation
# ---------------------------------------------------------------------------
print("[9] Lineage / Moller-key scoring ...")
lin = {}
with open(LIN, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        lin.setdefault(row["set"], []).append(row["probe"])

def mean_beta_full(probe_list, mat, p2i):
    cols = [p2i[p] for p in probe_list if p in p2i]
    if not cols:
        return np.full(mat.shape[0], np.nan)
    return mat[:, cols].mean(axis=1)

# Neou (full probe set for coverage)
sf1 = mean_beta_full(lin.get("SF1_lineage", []), Neou, neu_p2i)
tpit = mean_beta_full(lin.get("TPIT_lineage", []), Neou, neu_p2i)
pit1 = mean_beta_full(lin.get("PIT1_lineage", []), Neou, neu_p2i)
moller = mean_beta_full(lin.get("Moller_key", []), Neou, neu_p2i)
neu_lin = sf1 - (tpit + pit1) / 2.0
with open(os.path.join(OUT, "neou_lineage_scores.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample_id", "nfpt", "lineage", "sf1_score", "tpit_score", "pit1_score",
                "lineage_score", "moller_key_beta"])
    for i, s in enumerate(neu_samples):
        a = neu_anno.get(s, {})
        w.writerow([s, neu_nfpt[i], a.get("lineage", ""),
                    f"{sf1[i]:.4f}", f"{tpit[i]:.4f}", f"{pit1[i]:.4f}",
                    f"{neu_lin[i]:.4f}", f"{moller[i]:.4f}"])

# GSE207937 + GSE283928 lineage/Moller for cross-cohort NFPT comparison
sf1_7 = mean_beta_full(lin.get("SF1_lineage", []), X207, p207_i)
tpit_7 = mean_beta_full(lin.get("TPIT_lineage", []), X207, p207_i)
pit1_7 = mean_beta_full(lin.get("PIT1_lineage", []), X207, p207_i)
moller_7 = mean_beta_full(lin.get("Moller_key", []), X207, p207_i)
lin_7 = sf1_7 - (tpit_7 + pit1_7) / 2.0

sf1_2 = mean_beta_full(lin.get("SF1_lineage", []), X283, p283_i)
tpit_2 = mean_beta_full(lin.get("TPIT_lineage", []), X283, p283_i)
pit1_2 = mean_beta_full(lin.get("PIT1_lineage", []), X283, p283_i)
moller_2 = mean_beta_full(lin.get("Moller_key", []), X283, p283_i)
lin_2 = sf1_2 - (tpit_2 + pit1_2) / 2.0

# Neou NFPT vs Functioning
nf = neu_nfpt == 1; func = neu_nfpt == 0
neu_nfpt_lin = neu_lin[nf]; neu_func_lin = neu_lin[func]
neu_nfpt_mol = moller[nf]; neu_func_mol = moller[func]
# GSE207937 NS (NFPT) vs others
ns = gse207_nfpt == 1
g207_ns_lin = lin_7[ns]; g207_other_lin = lin_7[~ns]
g207_ns_mol = moller_7[ns]; g207_other_mol = moller_7[~ns]

print(f"    Neou lineage_score  NFPT={np.nanmean(neu_nfpt_lin):.3f}  Func={np.nanmean(neu_func_lin):.3f}")
print(f"    Neou Moller-key beta NFPT={np.nanmean(neu_nfpt_mol):.3f}  Func={np.nanmean(neu_func_mol):.3f}")
print(f"    GSE207937 lineage    NS={np.nanmean(g207_ns_lin):.3f}  other={np.nanmean(g207_other_lin):.3f}")
print(f"    GSE207937 Moller-key NS={np.nanmean(g207_ns_mol):.3f}  other={np.nanmean(g207_other_mol):.3f}")

# ---------------------------------------------------------------------------
# 10. Save summary
# ---------------------------------------------------------------------------
summary = {
    "neou_n_samples": int(Nc.shape[0]),
    "neou_n_nfpt": int(neu_nfpt.sum()),
    "neou_stage_counts": {k: neu_stage.count(k) for k in set(neu_stage)},
    "neou_lineage_counts": {k: neu_lineage.count(k) for k in set(neu_lineage)},
    "common_probes_3cohort": len(common),
    "variance_filtered_probes": int(mask.sum()),
    "pca_pc1": round(float(evr[0]), 4), "pca_pc2": round(float(evr[1]), 4),
    "pca_cum20": round(float(evr[:20].sum()), 4),
    "gse283928_tumor_clusters": {int(c): int((lab283 == c).sum()) for c in tumor_clusters},
    "neou_nfpt_vs_cluster": neu_ct,
    "threeway_cluster_counts": threeway,
    "neou_lineage_score_NFPT_mean": round(float(np.nanmean(neu_nfpt_lin)), 4),
    "neou_lineage_score_Func_mean": round(float(np.nanmean(neu_func_lin)), 4),
    "neou_lineage_score_perm_p": round(float(perm_p(neu_nfpt_lin, neu_func_lin)), 4),
    "neou_moller_beta_NFPT_mean": round(float(np.nanmean(neu_nfpt_mol)), 4),
    "neou_moller_beta_Func_mean": round(float(np.nanmean(neu_func_mol)), 4),
    "neou_moller_beta_perm_p": round(float(perm_p(neu_nfpt_mol, neu_func_mol)), 4),
    "gse207937_lineage_NS_mean": round(float(np.nanmean(g207_ns_lin)), 4),
    "gse207937_lineage_other_mean": round(float(np.nanmean(g207_other_lin)), 4),
    "gse207937_moller_NS_mean": round(float(np.nanmean(g207_ns_mol)), 4),
    "gse207937_moller_other_mean": round(float(np.nanmean(g207_other_mol)), 4),
}
with open(os.path.join(OUT, "neou_crosscohort_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
print(f"[done] total elapsed {time.time()-t0:.1f}s")
