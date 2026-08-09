#!/usr/bin/env python3
"""GSE207937 (Silva-Junior 2023) methylome stratification + NFPT biological validation.

Deep-dig goal: prove the MOLLER 2025 methylation subtype framework captures NFPT biology
by applying it to a SECOND, NFPT-labeled cohort (GSE207937 has subtype GH/ACTH/NF; NF=NFPT).
Also performs cross-cohort concordance vs GSE283928 (Belakhoua 2025).

Data notes:
- GSE207937 series matrix stores M-values (logit beta) -> convert to beta.
- GSE283928 stores beta directly. Both aligned to beta space for cross-cohort.
"""
import gzip, csv, os, json, time, re
import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from _paths import DATASET_PATH as BASE, DATASET_PATH_TEMP as TEMP
MAT_207937 = os.path.join(BASE, "GSE207937_series_matrix.txt.gz")
SAMPLES_207937 = os.path.join(TEMP, "samples", "GSE207937_samples.csv")
MAT_283928 = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
LINEAGE = os.path.join(TEMP, "results", "GSE283928", "lineage_marker_probes.csv")
CLUST_283928 = os.path.join(TEMP, "results", "GSE283928", "cluster_assignment_k5.csv")
OUT = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/results/GSE207937"
os.makedirs(OUT, exist_ok=True)

def clean(s): return s.strip().strip('"')

def m2beta(m): return 1.0 / (1.0 + np.power(2.0, -m))

t0 = time.time()

# ---------------------------------------------------------------------------
# 1. Parse GSE207937 series matrix -> beta matrix (probes x samples)
# ---------------------------------------------------------------------------
print("[1] Parsing GSE207937 matrix (M->beta) ...")
header = None
probes = []
rows = []
with gzip.open(MAT_207937, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!") or line.startswith("#") or line.strip() == "":
            continue
        parts = [clean(p) for p in line.rstrip("\n").split("\t")]
        if header is None:
            header = parts
            continue
        probes.append(parts[0])
        rows.append(parts[1:])
X = np.array(rows, dtype=np.float32)          # probes x samples, M-values
print(f"    raw shape {X.shape}, elapsed {time.time()-t0:.1f}s")
# M -> beta
Xb = m2beta(X).astype(np.float32)
del X
print(f"    beta range min={Xb.min():.3f} max={Xb.max():.3f} mean={Xb.mean():.3f}")
sample_ids = header[1:]
n_probe, n_samp = Xb.shape
print(f"    samples={n_samp}, probes={n_probe}")

# ---------------------------------------------------------------------------
# 2. Subtype labels (GH / ACTH / NF) from samples CSV
# ---------------------------------------------------------------------------
print("[2] Parsing subtype labels ...")
subtype_by_gsm = {}
with open(SAMPLES_207937, "r", encoding="utf-8", errors="ignore") as f:
    r = csv.DictReader(f)
    for row in r:
        gsm = clean(row.get("sample_geo_accession", ""))
        ch = row.get("sample_characteristics_ch1", "") or ""
        m = re.search(r"subtypes:\s*([A-Za-z\-]+)", ch)
        st = m.group(1) if m else "NA"
        if gsm:
            subtype_by_gsm[gsm] = st
# map sample_ids (GSMxxxx) -> subtype, build NFPT flag (NS = non-secreting = NFPT)
subtypes = []
nfpt = []
for s in sample_ids:
    st = subtype_by_gsm.get(s, "NA")
    subtypes.append(st)
    nfpt.append(1 if st == "NS" else 0)
print("    subtype counts:", {k: subtypes.count(k) for k in set(subtypes)})

# ---------------------------------------------------------------------------
# 3. Unsupervised 5-cluster reproducibility (variance filter -> PCA -> KMeans)
# ---------------------------------------------------------------------------
print("[3] Unsupervised 5-cluster ...")
var = Xb.var(axis=1)
k = max(2000, int(len(var) * 0.10))
thr = np.partition(var, -k)[-k]
mask = var >= thr
Xf = Xb[mask, :]
print(f"    variance-filtered probes kept: {int(mask.sum())}")
Xc = Xf - Xf.mean(axis=0, keepdims=True)
pca = PCA(n_components=20, random_state=42)
PC = pca.fit_transform(Xc.T)           # samples x 20
evr = pca.explained_variance_ratio_
print(f"    PC explained variance top5: {np.round(evr[:5],3)}  cum20={evr[:20].sum():.3f}")
km = KMeans(n_clusters=5, random_state=42, n_init=10)
labels = km.fit_predict(PC[:, :10])
uniq, cnts = np.unique(labels, return_counts=True)
sizes = dict(zip([int(u) for u in uniq], [int(c) for c in cnts]))
print(f"    cluster sizes: {sizes}")
sil = float(silhouette_score(PC[:, :10], labels))
print(f"    silhouette (k=5, top10PC): {sil:.3f}")

# between/total variance ratio as separation metric
centroids = km.cluster_centers_
# overall dispersion
def disp(P, lab):
    tot = 0.0; bet = 0.0; g = P.mean(axis=0)
    for c in np.unique(lab):
        m = lab == c; Pc = P[m]
        bet += len(Pc) * ((Pc.mean(axis=0) - g) ** 2).sum()
        tot += ((Pc - Pc.mean(axis=0)) ** 2).sum()
    return float(bet / tot) if tot > 0 else 0.0
sep_ratio = disp(PC[:, :10], labels)
print(f"    between/total variance ratio: {sep_ratio:.3f}")

# cluster x subtype contingency
print("[4] Cluster x subtype contingency ...")
subtypes_arr = np.array(subtypes)
contingency = {}
for c in np.unique(labels):
    m = labels == c
    sub = subtypes_arr[m]
    contingency[int(c)] = {k: int((sub == k).sum()) for k in ["GH", "ACTH", "NS", "NA"]}
print("   ", contingency)

# ---------------------------------------------------------------------------
# 4. Lineage score + Moller key-gene methylation (biological validation)
# ---------------------------------------------------------------------------
print("[5] Lineage / Moller-key methylation validation ...")
lin = {}
with open(LINEAGE, "r", encoding="utf-8", errors="ignore") as f:
    for row in csv.DictReader(f):
        lin.setdefault(row["set"], []).append(row["probe"])
# index probe->col
probe2col = {p: i for i, p in enumerate(probes)}
def mean_beta(probe_list):
    cols = [probe2col[p] for p in probe_list if p in probe2col]
    if not cols:
        return np.full(n_samp, np.nan)
    return Xb[cols, :].mean(axis=0)
sf1 = mean_beta(lin.get("SF1_lineage", []))
tpit = mean_beta(lin.get("TPIT_lineage", []))
pit1 = mean_beta(lin.get("PIT1_lineage", []))
moller = mean_beta(lin.get("Moller_key", []))
lineage_score = sf1 - (tpit + pit1) / 2.0
sf1_probes = lin.get("SF1_lineage", [])
sf1_cols = [p for p in sf1_probes if p in probe2col]
print(f"    SF1 probes={len(sf1_probes)} mapped={len(sf1_cols)}")
print(f"    lineage_score: NS(NFPT) mean={np.nanmean(lineage_score[np.array(nfpt)==1]):.3f}  secretory mean={np.nanmean(lineage_score[np.array(nfpt)==0]):.3f}")

# Moller key methylation by subtype
print("[6] Moller-key methylation by subtype ...")
moller_by_sub = {}
for st in ["GH", "ACTH", "NS"]:
    m = np.array(subtypes) == st
    moller_by_sub[st] = float(np.nanmean(moller[m])) if m.sum() > 0 else None
print("   ", moller_by_sub)

# ---------------------------------------------------------------------------
# 5. Cross-cohort concordance vs GSE283928 (Belakhoua)
# ---------------------------------------------------------------------------
print("[7] Cross-cohort: load GSE283928 + project GSE207937 onto its tumor-cluster centroids ...")
h2 = None; p2 = []; r2 = []
with gzip.open(MAT_283928, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!") or line.startswith("#") or line.strip() == "":
            continue
        parts = [clean(p) for p in line.rstrip("\n").split("\t")]
        if h2 is None:
            h2 = parts; continue
        p2.append(parts[0]); r2.append(parts[1:])
X2 = np.array(r2, dtype=np.float32)        # probes x samples (beta)
print(f"    GSE283928 shape {X2.shape}, beta range {X2.min():.3f}-{X2.max():.3f}")
s2 = h2[1:]
p2_idx = {p: i for i, p in enumerate(p2)}
# load cluster assignment, drop normal (cluster 1)
clu = {}
with open(CLUST_283928, "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        clu[row["sample_id"]] = int(row["cluster"])
# align columns of X2 to cluster labels
cl2 = np.array([clu.get(s, -1) for s in s2])
# common probes
set1 = set(probes); set2 = set(p2)
common = sorted(set1 & set2)
ci1 = [probe2col[p] for p in common]
ci2 = [p2_idx[p] for p in common]
X1c = Xb[ci1, :]      # common probes x 77 (GSE207937)
X2c = X2[ci2, :]      # common probes x 236 (GSE283928)
print(f"    common probes = {len(common)}")
# tumor-cluster centroids in common beta space
tumor_clusters = [c for c in np.unique(cl2) if c != 1 and c != -1]
centroids = {}
for c in tumor_clusters:
    m = cl2 == c
    centroids[int(c)] = X2c[:, m].mean(axis=1)
Cmat = np.array([centroids[c] for c in tumor_clusters])   # nclust x probes
# assign each GSE207937 sample to nearest centroid (euclidean), per-sample to limit memory
D = np.zeros((X1c.shape[1], len(tumor_clusters)), dtype=np.float32)
for i in range(X1c.shape[1]):
    D[i] = np.linalg.norm(Cmat - X1c[:, i:i+1].T, axis=1)
assigned = np.array(tumor_clusters)[np.argmin(D, axis=1)]
# also best correlation
cn = np.linalg.norm(Cmat, axis=1)
xn = np.linalg.norm(X1c, axis=0)
corr = (X1c.T @ Cmat.T) / (xn[:, None] * cn[None, :])
assigned_corr = np.array(tumor_clusters)[np.argmax(corr, axis=1)]
# cross-tab: GSE207937 subtype vs assigned GSE283928 tumor cluster
print("[8] Cross-cohort assignment (euclidean) by subtype ...")
xtab = {}
for st in ["GH", "ACTH", "NS"]:
    m = np.array(subtypes) == st
    vals, cc = np.unique(assigned[m], return_counts=True)
    xtab[st] = {int(v): int(c) for v, c in zip(vals, cc)}
print("   ", xtab)

# ---------------------------------------------------------------------------
# 6. Save outputs
# ---------------------------------------------------------------------------
print("[9] Saving outputs ...")
with open(os.path.join(OUT, "sample_labels.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample_id", "subtype", "nfpt", "cluster", "lineage_score", "moller_key_beta"])
    for i, s in enumerate(sample_ids):
        w.writerow([s, subtypes[i], nfpt[i], int(labels[i]),
                    f"{lineage_score[i]:.4f}", f"{moller[i]:.4f}"])
with open(os.path.join(OUT, "cluster_assignment_k5.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample_id", "cluster"])
    for s, c in zip(sample_ids, labels): w.writerow([s, int(c)])
with open(os.path.join(OUT, "cluster_subtype_contingency.json"), "w") as f:
    json.dump(contingency, f, indent=2)
with open(os.path.join(OUT, "moller_key_by_subtype.json"), "w") as f:
    json.dump(moller_by_sub, f, indent=2)
with open(os.path.join(OUT, "crosscohort_assignment.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["sample_id", "subtype", "assigned_cluster_eucl", "assigned_cluster_corr"])
    for i, s in enumerate(sample_ids):
        w.writerow([s, subtypes[i], int(assigned[i]), int(assigned_corr[i])])
summary = {
    "n_samples": int(n_samp), "n_probes": int(n_probe),
    "subtype_counts": {k: subtypes.count(k) for k in set(subtypes)},
    "cluster_sizes": sizes,
    "silhouette_k5": round(sil, 4),
    "between_total_var_ratio": round(sep_ratio, 4),
    "pc_top5_explained": [round(float(x), 3) for x in evr[:5]],
    "lineage_score_NS_mean": round(float(np.nanmean(lineage_score[np.array(nfpt) == 1])), 4),
    "lineage_score_secretory_mean": round(float(np.nanmean(lineage_score[np.array(nfpt) == 0])), 4),
    "moller_key_beta_by_subtype": moller_by_sub,
    "crosscohort_n_common_probes": len(common),
    "crosscohort_tumor_clusters": [int(c) for c in tumor_clusters],
    "crosscohort_subtype_vs_cluster": xtab,
}
with open(os.path.join(OUT, "analysis_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
print(f"[done] total elapsed {time.time()-t0:.1f}s")
