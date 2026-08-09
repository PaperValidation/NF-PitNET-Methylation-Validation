#!/usr/bin/env python3
"""Neou 2020 (E-MTAB-7762) INTERNAL independent clustering (Option C).

Goal: recapitulate the methylation taxonomy FROM SCRATCH on Neou alone (86 EPIC
samples), without using the GSE283928 centroid projection. Then validate that the
internal structure agrees with (a) the cross-cohort projection (neu_e onto GSE283928
tumor clusters) and (b) Møller-2025 key-gene hypermethylation / lineage scores.

Design (mirrors neou_crosscohort.py for comparability):
  - load Neou + GSE283928(tumor) + GSE207937 in common-probe beta space
  - REFERENCE space: combined variance filter (top 10%) + GSE283928 centroids (reused
    for cross-cohort mapping of each Neou sample)
  - INTERNAL space: Neou-own variance filter (top 10%) + PCA(20, seed=42) -> KMeans k=2..6
  - pick k by silhouette; characterize clusters by NFPT fraction, lineage score,
    Møller-key beta, and cross-tab vs cross-cohort assigned cluster
Outputs to results/Neou/:
  neou_internal_clusters.csv, neou_internal_silhouette.json,
  neou_internal_cluster_summary.json, neou_internal_cluster_figure.png
"""
import gzip, csv, json, time, re, os
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from _paths import DATASET_PATH as BASE, DATASET_PATH_TEMP as TEMP
NE = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
MAT_283 = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
MAT_207 = os.path.join(BASE, "GSE207937_series_matrix.txt.gz")
SAMP_207 = os.path.join(TEMP, "samples", "GSE207937_samples.csv")
CLUST_283 = os.path.join(TEMP, "results", "GSE283928", "cluster_assignment_k5.csv")
LIN = os.path.join(TEMP, "results", "GSE283928", "lineage_marker_probes.csv")
NEU_ANNO = os.path.join(BASE, "E-MTAB-7762", "beta_out", "neou_sample_annotation.csv")
NEU_BIN = os.path.join(BASE, "E-MTAB-7762", "beta_out", "neou_beta_full.bin")
NEU_PROBES = os.path.join(BASE, "E-MTAB-7762", "beta_out", "neou_beta_probes.txt")
NEU_SAMPLES = os.path.join(BASE, "E-MTAB-7762", "beta_out", "neou_beta_samples.txt")
OUT = os.path.join(NE, "results", "Neou")
os.makedirs(OUT, exist_ok=True)

def clean(s): return s.strip().strip('"')
def m2beta(m): return 1.0 / (1.0 + np.power(2.0, -m))
def impute(M):
    col_mean = np.nanmean(M, axis=0)
    col_mean = np.where(np.isnan(col_mean), 0.5, col_mean)
    idx = np.where(np.isnan(M))
    if len(idx[0]):
        M = M.copy(); M[idx] = np.take(col_mean, idx[1])
    return M
def perm_p(a, b, n=20000, seed=0):
    rng = np.random.default_rng(seed)
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) == 0 or len(b) == 0: return float("nan")
    obs = a.mean() - b.mean(); pool = np.concatenate([a, b]); cnt = 0
    for _ in range(n):
        rng.shuffle(pool); d = pool[:len(a)].mean() - pool[len(a):].mean()
        if abs(d) >= abs(obs): cnt += 1
    return (cnt + 1) / (n + 1)
def topk_var_mask(M, frac=0.10):
    var = M.var(axis=0); k = max(1000, int(M.shape[1] * frac))
    thr = np.partition(var, -k)[-k]; return var >= thr

t0 = time.time()
print("[1] Load Neou beta ...")
with open(NEU_BIN, "rb") as f:
    nrow = int(np.fromfile(f, dtype=np.int32, count=1)[0])
    ncol = int(np.fromfile(f, dtype=np.int32, count=1)[0])
    data = np.fromfile(f, dtype=np.float64, count=nrow * ncol)
Bmat = data.reshape(nrow, ncol, order="F")
Neou = np.ascontiguousarray(Bmat.T).astype(np.float32)   # 86 x 866238
with open(NEU_PROBES, encoding="utf-8") as f: neu_probes = [l.rstrip("\n") for l in f]
with open(NEU_SAMPLES, encoding="utf-8") as f: neu_samples = [l.rstrip("\n") for l in f]
neu_p2i = {p: i for i, p in enumerate(neu_probes)}
Neou = impute(Neou)
print(f"    Neou {Neou.shape}; elapsed {time.time()-t0:.1f}s")

print("[2] Neou annotation ...")
neu_anno = {}
with open(NEU_ANNO, encoding="utf-8") as f:
    for row in csv.DictReader(f): neu_anno[row["individual"]] = row
neu_nfpt = np.array([int(neu_anno.get(s, {}).get("nfpt", 0)) for s in neu_samples])
neu_stage = [neu_anno.get(s, {}).get("disease_staging", "NA") for s in neu_samples]
neu_lineage = [neu_anno.get(s, {}).get("lineage", "NA") for s in neu_samples]
print(f"    NFPT={int(neu_nfpt.sum())}/{len(neu_samples)}")

print("[3] Load GSE283928 tumor + GSE207937 (for reference-space projection) ...")
with gzip.open(MAT_283, "rt", encoding="utf-8", errors="ignore") as f:
    hdr = f.readline().rstrip("\n").split("\t")
beta_idx = [i for i, h in enumerate(hdr) if h.endswith("_Beta")]
beta_cols = [hdr[i] for i in beta_idx]
rows = []; probes283 = []
with gzip.open(MAT_283, "rt", encoding="utf-8", errors="ignore") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t"); probes283.append(p[0]); rows.append([p[i] for i in beta_idx])
X283 = impute(np.array(rows, dtype=np.float32).T); p283 = np.array(probes283); p283_i = {p: j for j, p in enumerate(p283)}
del rows
clu = {}
with open(CLUST_283, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["sample_id"].endswith("_Beta"): clu[row["sample_id"]] = int(row["cluster"])
lab283 = np.array([clu.get(c, -1) for c in beta_cols])
tumor_clusters = sorted(set(lab283.tolist()) - {-1})

h2 = None; p207 = []; r207 = []
with gzip.open(MAT_207, "rt", encoding="utf-8", errors="ignore") as f:
    for line in f:
        if line.startswith("!") or line.startswith("#") or not line.strip(): continue
        parts = [clean(x) for x in line.rstrip("\n").split("\t")]
        if h2 is None: h2 = parts; continue
        p207.append(parts[0]); r207.append(parts[1:])
X207 = impute(m2beta(np.array(r207, dtype=np.float32)).T); p207 = np.array(p207); p207_i = {p: j for j, p in enumerate(p207)}
gse207_samples = h2[1:]
sub_by_gsm = {}
with open(SAMP_207, encoding="utf-8", errors="ignore") as f:
    for row in csv.DictReader(f):
        gsm = clean(row.get("sample_geo_accession", "")); ch = row.get("sample_characteristics_ch1", "") or ""
        m = re.search(r"subtypes:\s*([A-Za-z\-]+)", ch); sub_by_gsm[gsm] = m.group(1) if m else "NA"
gse207_sub = [sub_by_gsm.get(s, "NA") for s in gse207_samples]
gse207_nfpt = np.array([1 if s == "NS" else 0 for s in gse207_sub])

print("[4] Common probes + REFERENCE space (combined variance filter) ...")
set_n = set(neu_probes); set_2 = set(p283.tolist()); set_7 = set(p207.tolist())
common = sorted(set_n & set_2 & set_7)
ci_n = [neu_p2i[p] for p in common]; ci_2 = [p283_i[p] for p in common]; ci_7 = [p207_i[p] for p in common]
Nc = Neou[:, ci_n]; X2c = X283[:, ci_2]; X7c = X207[:, ci_7]
combined = np.vstack([Nc, X2c, X7c]).astype(np.float32)
ref_mask = topk_var_mask(combined, 0.10)
Xref = combined[:, ref_mask]                       # reference space (all 3 cohorts, for record)
# GSE283928 centroids in the UNFILTERED common-probe beta space (matches neou_crosscohort.py)
Cmat = np.zeros((len(tumor_clusters), X2c.shape[1]), dtype=np.float32)
for j, c in enumerate(tumor_clusters): Cmat[j] = X2c[lab283 == c].mean(axis=0)
cn = np.linalg.norm(Cmat, axis=1)
def assign(Xq):
    D = np.zeros((Xq.shape[0], len(tumor_clusters)), dtype=np.float32)
    for i in range(Xq.shape[0]): D[i] = np.linalg.norm(Cmat - Xq[i:i+1], axis=1)
    return np.array(tumor_clusters)[np.argmin(D, axis=1)], D
neu_e, neu_D = assign(Nc)                          # cross-cohort projection (reference)
print(f"    common={len(common)}; ref_mask kept {int(ref_mask.sum())}; Neou cross-cohort assign:",
      {int(c): int((neu_e == c).sum()) for c in tumor_clusters})

# lineage / moller on full Neou probe set
lin = {}
with open(LIN, encoding="utf-8") as f:
    for row in csv.DictReader(f): lin.setdefault(row["set"], []).append(row["probe"])
def mb(probe_list, mat, p2i):
    cols = [p2i[p] for p in probe_list if p in p2i]
    return np.full(mat.shape[0], np.nan) if not cols else mat[:, cols].mean(axis=1)
sf1 = mb(lin.get("SF1_lineage", []), Neou, neu_p2i); tpit = mb(lin.get("TPIT_lineage", []), Neou, neu_p2i)
pit1 = mb(lin.get("PIT1_lineage", []), Neou, neu_p2i); moller = mb(lin.get("Moller_key", []), Neou, neu_p2i)
neu_lin = sf1 - (tpit + pit1) / 2.0

print("[5] INTERNAL space: Neou-own variance filter + PCA(20) ...")
int_mask = topk_var_mask(Neou, 0.10)               # Neou-INTERNAL top-10% variable probes
Xint = Neou[:, int_mask].astype(np.float32)
pca = PCA(n_components=20, random_state=42); PC = pca.fit_transform(Xint)
evr = pca.explained_variance_ratio_
print(f"    internal kept {int(int_mask.sum())} probes; PC1/PC2={evr[0]:.3f}/{evr[1]:.3f}; cum20={evr[:20].sum():.3f}")

print("[6] KMeans k=2..6 + silhouette ...")
sil = {}
best_k, best_s = None, -1
for k in range(2, 7):
    km = KMeans(n_clusters=k, random_state=42, n_init=10).fit(PC)
    s = silhouette_score(PC, km.labels_)
    sil[k] = round(float(s), 4)
    if s > best_s: best_s, best_k = s, k
    print(f"    k={k} silhouette={s:.4f}")
print(f"    best_k={best_k} (silhouette {best_s:.4f})")

# final clustering at best_k
km = KMeans(n_clusters=best_k, random_state=42, n_init=10).fit(PC)
internal = km.labels_
# order clusters by mean Møller-key beta (high = high-risk-like) for stable labeling
order = np.argsort([moller[internal == c].mean() for c in range(best_k)])[::-1]
remap = {old: new for new, old in enumerate(order)}
internal = np.array([remap[l] for l in internal])   # cluster 0 = highest Møller key beta

print("[7] Characterize internal clusters ...")
char = {}
for c in range(best_k):
    m = internal == c
    xproj = {int(cl): int((neu_e[m] == cl).sum()) for cl in tumor_clusters}
    char[c] = {
        "n": int(m.sum()),
        "nfpt": int(neu_nfpt[m].sum()),
        "nfpt_frac": round(float(neu_nfpt[m].mean()), 3),
        "moller_key_beta": round(float(np.nanmean(moller[m])), 4),
        "lineage_score": round(float(np.nanmean(neu_lin[m])), 4),
        "cross_cohort_cluster_counts": xproj,
        "stages": {k: int((np.array(neu_stage)[m] == k).sum()) for k in set(np.array(neu_stage)[m])},
        "lineages": {k: int((np.array(neu_lineage)[m] == k).sum()) for k in set(np.array(neu_lineage)[m])},
    }
    print(f"    cluster {c}: n={char[c]['n']} NFPT={char[c]['nfpt']}({char[c]['nfpt_frac']}) "
          f"Mollerβ={char[c]['moller_key_beta']} lin={char[c]['lineage_score']} xproj={xproj}")

# cross-tab internal x nfpt (does internal clustering recapitulate NFPT/Functioning split?)
print("[8] Internal-cluster x NFPT contingency ...")
ct = {}
for st, fl in [("NFPT", 1), ("Functioning", 0)]:
    mm = neu_nfpt == fl
    vals, cc = np.unique(internal[mm], return_counts=True)
    ct[st] = {int(v): int(n) for v, n in zip(vals, cc)}
print("   ", ct)

# association test: internal cluster (2-group if best_k==2) vs NFPT via permutation
if best_k == 2:
    grp0 = neu_nfpt[internal == 0]; grp1 = neu_nfpt[internal == 1]
    assoc_p = perm_p(grp0, grp1)
    print(f"    NFPT-rate cluster0={grp0.mean():.3f} cluster1={grp1.mean():.3f} perm_p={assoc_p:.4g}")

# save outputs
with open(os.path.join(OUT, "neou_internal_clusters.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["sample_id", "disease_staging", "nfpt", "lineage", "internal_cluster",
                "cross_cohort_cluster_eucl", "lineage_score", "moller_key_beta"])
    for i, s in enumerate(neu_samples):
        a = neu_anno.get(s, {})
        w.writerow([s, a.get("disease_staging", ""), neu_nfpt[i], a.get("lineage", ""),
                    int(internal[i]), int(neu_e[i]), f"{neu_lin[i]:.4f}", f"{moller[i]:.4f}"])
with open(os.path.join(OUT, "neou_internal_silhouette.json"), "w") as f:
    json.dump({"silhouette_by_k": sil, "best_k": int(best_k), "best_silhouette": round(float(best_s), 4)}, f, indent=2)
summary = {
    "neou_n_samples": int(Neou.shape[0]),
    "neou_n_nfpt": int(neu_nfpt.sum()),
    "internal_n_probes": int(int_mask.sum()),
    "pca_pc1": round(float(evr[0]), 4), "pca_pc2": round(float(evr[1]), 4), "pca_cum20": round(float(evr[:20].sum()), 4),
    "best_k": int(best_k), "silhouette_by_k": sil,
    "cluster_characterization": char,
    "internal_x_nfpt_contingency": ct,
    "assoc_perm_p": round(float(assoc_p), 6) if best_k == 2 else None,
    "cross_cohort_assignment_counts": {int(c): int((neu_e == c).sum()) for c in tumor_clusters},
}
with open(os.path.join(OUT, "neou_internal_cluster_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))

# ---------- figure ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 2, figsize=(11, 9))
# (0,0) PCA scatter colored by internal cluster, NFPT shape
shapes = ["o" if neu_nfpt[i] == 1 else "s" for i in range(len(neu_samples))]
cmap = plt.cm.tab10
for c in range(best_k):
    m = internal == c
    ax[0, 0].scatter(PC[m, 0], PC[m, 1], s=46, alpha=0.8,
                     label=f"int{c} (n={int(m.sum())})", color=cmap(c % 10))
for i in range(PC.shape[0]):
    pass
ax[0, 0].set_xlabel(f"PC1 ({evr[0]*100:.1f}%)"); ax[0, 0].set_ylabel(f"PC2 ({evr[1]*100:.1f}%)")
ax[0, 0].set_title("Neou internal PCA — KMeans (from scratch)"); ax[0, 0].legend(fontsize=7)

# (0,1) silhouette vs k
ks = sorted(sil.keys()); ax[0, 1].plot(ks, [sil[k] for k in ks], "o-", color="steelblue")
ax[0, 1].axvline(best_k, ls="--", color="red", alpha=0.6); ax[0, 1].set_xlabel("k"); ax[0, 1].set_ylabel("silhouette")
ax[0, 1].set_title(f"Best k = {best_k}")

# (1,0) NFPT fraction + Møller key beta per internal cluster
cs = list(range(best_k))
ax[1, 0].bar(cs, [char[c]["nfpt_frac"] for c in cs], color=[cmap(c % 10) for c in cs], alpha=0.85)
ax[1, 0].set_xlabel("internal cluster"); ax[1, 0].set_ylabel("NFPT fraction", color="darkgreen")
ax2 = ax[1, 0].twinx()
ax2.plot(cs, [char[c]["moller_key_beta"] for c in cs], "D-", color="darkred")
ax2.set_ylabel("Møller-key β", color="darkred"); ax[1, 0].set_title("NFPT fraction & Møller-key β by cluster")

# (1,1) stacked cross-cohort projection per internal cluster
bottom = np.zeros(best_k)
for j, cl in enumerate(tumor_clusters):
    vals = [char[c]["cross_cohort_cluster_counts"].get(int(cl), 0) for c in cs]
    ax[1, 1].bar(cs, vals, bottom=bottom, label=f"->c{cl}")
    bottom += np.array(vals)
ax[1, 1].set_xlabel("internal cluster"); ax[1, 1].set_ylabel("Neou samples"); ax[1, 1].set_title("Internal cluster -> cross-cohort projection")
ax[1, 1].legend(fontsize=7)
plt.tight_layout()
fig.savefig(os.path.join(OUT, "neou_internal_cluster_figure.png"), dpi=150)
print(f"[done] elapsed {time.time()-t0:.1f}s; figure saved")
