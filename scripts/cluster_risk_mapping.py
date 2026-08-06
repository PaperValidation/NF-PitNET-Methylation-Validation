#!/usr/bin/env python3
"""Audit fix #1: GSE283928 cluster -> Moller risk hierarchy mapping.

QUESTION: NFPT maps to GSE283928 cluster 0. Is cluster 0 the HIGH-risk
cluster (Moller k3/k4/k5, HR 2.57-3.63) or the LOW-risk cluster (k1/k2)?

APPROACH:
  1. Per-cluster Moller-key gene mean beta (168 probes) — high beta = high risk
  2. Per-cluster lineage scores (SF1 / TPIT / PIT1) — identifies the TPIT/PIT1 cluster
  3. Per-cluster global mean beta (control for global hypermethylation)
  4. Pairwise permutation p-values between clusters
  5. Leave-one-out centroid assignment on GSE283928 itself (centroid bias check)
  6. Permutation null for "fraction assigned to c0" under random label shuffle
"""
import gzip, csv, json, time, os
import numpy as np

BASE = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
NE = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
MAT_283 = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
CLUST_283 = os.path.join(BASE, "results", "GSE283928", "cluster_assignment_k5.csv")
LIN = os.path.join(BASE, "results", "GSE283928", "lineage_marker_probes.csv")
OUT = os.path.join(NE, "results", "Neou")

t0 = time.time()

# -----------------------------------------------------------------------
# 1. Load GSE283928 tumor beta matrix (Beta cols only)
# -----------------------------------------------------------------------
print("[1] Loading GSE283928 (Beta cols only) ...")
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
X = np.array(rows, dtype=np.float32).T   # samples(118) x probes
del rows
p283 = np.array(probes283)
p283_i = {p: j for j, p in enumerate(p283)}
print(f"    shape {X.shape}; elapsed {time.time()-t0:.1f}s")

# impute NaN
nan_count = np.isnan(X).sum()
print(f"    NaN before impute: {nan_count} ({100*nan_count/X.size:.2f}%)")
col_mean = np.nanmean(X, axis=0)
col_mean = np.where(np.isnan(col_mean), 0.5, col_mean)
idx = np.where(np.isnan(X))
if len(idx[0]):
    X[idx] = np.take(col_mean, idx[1])

# -----------------------------------------------------------------------
# 2. Cluster labels (tumor only)
# -----------------------------------------------------------------------
print("[2] Cluster labels ...")
clu = {}
with open(CLUST_283, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["sample_id"].endswith("_Beta"):
            clu[row["sample_id"]] = int(row["cluster"])
lab = np.array([clu.get(c, -1) for c in beta_cols])
tumor_clusters = sorted(set(lab.tolist()) - {-1})
print(f"    tumor clusters {tumor_clusters}; sizes",
      {int(c): int((lab == c).sum()) for c in tumor_clusters})

# -----------------------------------------------------------------------
# 3. Load probe sets
# -----------------------------------------------------------------------
print("[3] Loading probe sets ...")
lin = {}
with open(LIN, encoding="utf-8") as f:
    for row in csv.DictReader(f):
        lin.setdefault(row["set"], []).append(row["probe"])
moller_probes = lin.get("Moller_key", [])
sf1_probes = lin.get("SF1_lineage", [])
tpit_probes = lin.get("TPIT_lineage", [])
pit1_probes = lin.get("PIT1_lineage", [])
print(f"    Moller_key={len(moller_probes)} SF1={len(sf1_probes)} "
      f"TPIT={len(tpit_probes)} PIT1={len(pit1_probes)}")

def mean_beta(probe_list, mat, p2i):
    cols = [p2i[p] for p in probe_list if p in p2i]
    if not cols:
        return np.full(mat.shape[0], np.nan)
    return mat[:, cols].mean(axis=1)

# per-sample scores
moller = mean_beta(moller_probes, X, p283_i)
sf1 = mean_beta(sf1_probes, X, p283_i)
tpit = mean_beta(tpit_probes, X, p283_i)
pit1 = mean_beta(pit1_probes, X, p283_i)
lineage_score = sf1 - (tpit + pit1) / 2.0
global_beta = X.mean(axis=1)

# -----------------------------------------------------------------------
# 4. Per-cluster statistics
# -----------------------------------------------------------------------
print("\n[4] Per-cluster statistics:")
print(f"    {'Cluster':>8} {'N':>5} {'Moller_beta':>12} {'SF1':>8} {'TPIT':>8} "
      f"{'PIT1':>8} {'LinScore':>9} {'GlobalBeta':>11}")
print(f"    {'-'*75}")
cluster_stats = {}
for c in tumor_clusters:
    m = lab == c
    n = int(m.sum())
    mb = float(moller[m].mean())
    s = float(sf1[m].mean())
    t = float(tpit[m].mean())
    p = float(pit1[m].mean())
    ls = float(lineage_score[m].mean())
    gb = float(global_beta[m].mean())
    # SEM
    mb_sem = float(moller[m].std() / np.sqrt(n)) if n > 1 else 0.0
    cluster_stats[int(c)] = {
        "n": n, "moller_beta_mean": round(mb, 4), "moller_beta_sem": round(mb_sem, 4),
        "sf1_mean": round(s, 4), "tpit_mean": round(t, 4), "pit1_mean": round(p, 4),
        "lineage_score": round(ls, 4), "global_beta_mean": round(gb, 4)
    }
    print(f"    {'c'+str(c):>8} {n:>5} {mb:>12.4f} {s:>8.4f} {t:>8.4f} "
          f"{p:>8.4f} {ls:>9.4f} {gb:>11.4f}")

# overall
print(f"    {'ALL':>8} {len(lab):>5} {moller.mean():>12.4f} {sf1.mean():>8.4f} "
      f"{tpit.mean():>8.4f} {pit1.mean():>8.4f} {lineage_score.mean():>9.4f} "
      f"{global_beta.mean():>11.4f}")

# -----------------------------------------------------------------------
# 5. Pairwise permutation p-values (Moller key beta)
# -----------------------------------------------------------------------
print("\n[5] Pairwise permutation p-values (Moller key beta):")
def perm_p(a, b, n=10000, seed=0):
    rng = np.random.default_rng(seed)
    a = np.asarray(a, float); b = np.asarray(b, float)
    obs = a.mean() - b.mean()
    pool = np.concatenate([a, b])
    cnt = 0
    for _ in range(n):
        rng.shuffle(pool)
        d = pool[:len(a)].mean() - pool[len(a):].mean()
        if abs(d) >= abs(obs):
            cnt += 1
    return (cnt + 1) / (n + 1)

pairwise = {}
for i, ci in enumerate(tumor_clusters):
    for cj in tumor_clusters[i+1:]:
        mi = lab == ci; mj = lab == cj
        p = perm_p(moller[mi], moller[mj])
        diff = float(moller[mi].mean() - moller[mj].mean())
        pairwise[f"c{ci}_vs_c{cj}"] = {"diff": round(diff, 4), "perm_p": round(p, 5)}
        sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
        print(f"    c{ci} vs c{cj}: diff={diff:+.4f}  p={p:.5f}  {sig}")

# -----------------------------------------------------------------------
# 6. Leave-one-out centroid assignment (centroid bias check)
# -----------------------------------------------------------------------
print("\n[6] Leave-one-out centroid assignment (bias check):")
# Use common-probe space (all probes here since single-cohort)
loo_assigned = np.zeros(len(lab), dtype=int)
loo_dist_to_c0 = np.zeros(len(lab))
for i in range(len(lab)):
    mask = np.ones(len(lab), dtype=bool); mask[i] = False
    centroids = np.array([X[mask & (lab == c)].mean(axis=0) for c in tumor_clusters])
    d = np.linalg.norm(centroids - X[i:i+1], axis=1)
    loo_assigned[i] = tumor_clusters[np.argmin(d)]
    loo_dist_to_c0[i] = d[tumor_clusters.index(0)] if 0 in tumor_clusters else np.nan

loo_correct = np.mean(loo_assigned == lab)
loo_c0_frac = np.mean(loo_assigned == 0)
print(f"    LOO accuracy: {loo_correct:.3f} ({int((loo_assigned==lab).sum())}/{len(lab)})")
print(f"    LOO fraction assigned to c0: {loo_c0_frac:.3f}")
loo_ct = {}
for c in tumor_clusters:
    m = lab == c
    vals, cc = np.unique(loo_assigned[m], return_counts=True)
    loo_ct[int(c)] = {int(v): int(cnt) for v, cnt in zip(vals, cc)}
print(f"    LOO confusion: {loo_ct}")

# -----------------------------------------------------------------------
# 7. Permutation null: shuffle NFPT labels, fraction->c0
# -----------------------------------------------------------------------
print("\n[7] Permutation null for NFPT->c0 enrichment:")
# Load Neou assignment results
neu_assign = {}
with open(os.path.join(OUT, "neou_crosscohort_assignment.csv"), encoding="utf-8") as f:
    for row in csv.DictReader(f):
        neu_assign[row["sample_id"]] = int(row["assigned_cluster_eucl"])
neu_nfpt = {}
with open(os.path.join(OUT, "neou_crosscohort_assignment.csv"), encoding="utf-8") as f:
    for row in csv.DictReader(f):
        neu_nfpt[row["sample_id"]] = int(row["nfpt"])

nfpt_ids = [s for s, v in neu_nfpt.items() if v == 1]
n_nfpt = len(nfpt_ids)
all_assign = np.array([neu_assign[s] for s in neu_assign])
nfpt_assign = np.array([neu_assign[s] for s in nfpt_ids])
observed_c0_frac = np.mean(nfpt_assign == 0)
print(f"    Neou NFPT n={n_nfpt}; observed fraction->c0 = {observed_c0_frac:.3f} "
      f"({int((nfpt_assign==0).sum())}/{n_nfpt})")

# null: randomly pick n_nfpt samples from all 86, compute fraction->c0
rng = np.random.default_rng(42)
n_perm = 10000
null_c0_frac = np.zeros(n_perm)
for i in range(n_perm):
    idx = rng.choice(len(all_assign), size=n_nfpt, replace=False)
    null_c0_frac[i] = np.mean(all_assign[idx] == 0)
perm_p_nfpt = np.mean(null_c0_frac >= observed_c0_frac)
print(f"    Null (random {n_nfpt} from 86): mean={null_c0_frac.mean():.3f}, "
      f"std={null_c0_frac.std():.3f}")
print(f"    Permutation p (one-sided): {perm_p_nfpt:.5f}")
ci_lo, ci_hi = np.percentile(null_c0_frac, [2.5, 97.5])
print(f"    95% null CI: [{ci_lo:.3f}, {ci_hi:.3f}]")

# -----------------------------------------------------------------------
# 8. Save results
# -----------------------------------------------------------------------
results = {
    "per_cluster_stats": cluster_stats,
    "pairwise_moller_perm_p": pairwise,
    "loo_accuracy": round(float(loo_correct), 4),
    "loo_fraction_to_c0": round(float(loo_c0_frac), 4),
    "loo_confusion": loo_ct,
    "neou_nfpt_observed_c0_frac": round(float(observed_c0_frac), 4),
    "neou_nfpt_null_c0_frac_mean": round(float(null_c0_frac.mean()), 4),
    "neou_nfpt_null_c0_frac_std": round(float(null_c0_frac.std()), 4),
    "neou_nfpt_null_95ci": [round(float(ci_lo), 4), round(float(ci_hi), 4)],
    "neou_nfpt_perm_p": round(float(perm_p_nfpt), 5),
    "interpretation": ""
}

# Auto-interpretation
c0_moller = cluster_stats[0]["moller_beta_mean"]
max_moller_cluster = max(cluster_stats, key=lambda k: cluster_stats[k]["moller_beta_mean"])
c0_is_highest = (max_moller_cluster == 0)
c0_lineage = cluster_stats[0]["lineage_score"]
min_lin_cluster = min(cluster_stats, key=lambda k: cluster_stats[k]["lineage_score"])

interp = (
    f"Cluster 0 Moller-key beta = {c0_moller:.4f}. "
    f"Highest cluster = c{max_moller_cluster} ({cluster_stats[max_moller_cluster]['moller_beta_mean']:.4f}). "
    f"Cluster 0 IS {'the highest' if c0_is_highest else 'NOT the highest'} -> "
    f"{'NFPT maps to HIGH-risk cluster (narrative SUPPORTED)' if c0_is_highest else 'NFPT maps to LOW/MID cluster (narrative needs revision)'}. "
    f"Lowest lineage score cluster = c{min_lin_cluster} (lineage_score={cluster_stats[min_lin_cluster]['lineage_score']:.4f}) "
    f"-> this is the TPIT/PIT1 cluster (Moller k4 equivalent). "
    f"LOO centroid accuracy = {loo_correct:.3f} (bias check). "
    f"NFPT->c0 permutation p = {perm_p_nfpt:.5f} "
    f"({'significant' if perm_p_nfpt < 0.05 else 'NOT significant'})."
)
results["interpretation"] = interp
print(f"\n[8] Interpretation:\n    {interp}")

with open(os.path.join(OUT, "cluster_risk_mapping.json"), "w") as f:
    json.dump(results, f, indent=2)
print(f"\n[done] saved cluster_risk_mapping.json; elapsed {time.time()-t0:.1f}s")
