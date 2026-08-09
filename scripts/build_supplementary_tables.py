#!/usr/bin/env python3
# ============================================================================
# Build Supplementary Tables 1, 2, 3 for the NF-PitNET manuscript.
#   ST1 : Joint PCA cohort centroids + explained variance (PC1-20)
#   ST2 : Neou internal-cluster characterization (de-novo k=3 taxonomy)
#   ST3 : Normalization-sensitivity (preprocessRaw / noob / funnorm / quantile)
# Outputs: results/Neou/supplementary_table1_pca.csv
#          results/Neou/supplementary_table2_internal_clusters.csv
#          results/Neou/supplementary_table3_normalization.csv
# ============================================================================
import csv, json, os
import numpy as np

NE = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
OUT = os.path.join(NE, "results", "Neou")
from _paths import DATASET_PATH as G
# ST1 : Joint PCA cohort centroids + explained variance
# ---------------------------------------------------------------------------
pc = list(csv.DictReader(open(os.path.join(OUT, "joint_pca_coords.csv"))))
cohorts = ["Neou", "GSE283928", "GSE207937"]
pcs = [f"PC{i+1}" for i in range(20)]
M = {c: np.array([[float(r[p]) for p in pcs] for r in pc if r["cohort"] == c]) for c in cohorts}

# per-cohort centroid in PC space
cent_df = {"cohort": cohorts}
for p in pcs:
    cent_df[p] = [round(float(M[c][:, int(p[2:])-1].mean()), 4) for c in cohorts]
with open(os.path.join(OUT, "supplementary_table1_pca.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["table"] + ["component"] + ["Neou", "GSE283928", "GSE207937"])
    # part A: cohort centroids
    w.writerow(["A_cohort_centroid"] + ["(mean PC coord)"] + ["", "", ""])
    for i, p in enumerate(pcs):
        w.writerow(["A_cohort_centroid", p] + [cent_df[p][0], cent_df[p][1], cent_df[p][2]])
    # part B: explained variance (recovered exactly from PC score variance)
    evr = json.load(open(os.path.join(OUT, "neou_crosscohort_summary.json")))
    cum20 = float(evr["pca_cum20"])
    # TRUE eigenvalue_i = variance of the POOLED PCA scores (PCA was fit on the
    # combined 281-sample matrix, so var over all rows == explained_variance_[i]).
    pooled = np.array([[float(r[p]) for p in pcs] for r in pc])  # 281 x 20
    eig = pooled.var(axis=0, ddof=1)                            # == eigenvalues
    S = eig.sum()
    ratio = eig * cum20 / S                                     # == explained_variance_ratio_
    w.writerow(["B_explained_variance"] + ["(fraction)"] + ["", "", ""])
    cum = 0.0
    for i in range(20):
        cum += ratio[i]
        w.writerow(["B_explained_variance", pcs[i]] + [round(float(ratio[i]), 4), round(float(cum), 4), ""])
print("[ST1] PC1 ratio=%.4f (expect ~0.2923), PC2=%.4f (expect ~0.0967), cum20=%.4f (expect 0.6344)"
      % (ratio[0], ratio[1], ratio[:20].sum()))

# GSE283928 per-cluster centroids in PCA space (enrichment: reference taxonomy separation)
cl = {r["sample_id"]: int(r["cluster"]) for r in csv.DictReader(open(os.path.join(G, "results", "GSE283928", "cluster_assignment_k5.csv")))}
gse = [r for r in pc if r["cohort"] == "GSE283928"]
clusters = [0, 2, 3, 4]
clrows = []
for cc in clusters:
    sub = np.array([[float(r[p]) for p in pcs] for r in gse if cl.get(r["sample_id"]) == cc])
    clrows.append((cc, sub.shape[0], sub.mean(axis=0)))
with open(os.path.join(OUT, "supplementary_table1_pca_cluster_centroids.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gse283928_cluster", "n"] + pcs)
    for cc, n, mu in clrows:
        w.writerow([f"c{cc}", n] + [f"{x:.4f}" for x in mu])
print("[ST1] GSE283928 cluster centroids: c0 n=%d, c2 n=%d, c3 n=%d, c4 n=%d" % tuple(r[1] for r in clrows))

# ---------------------------------------------------------------------------
# ST2 : Neou internal-cluster characterization
# ---------------------------------------------------------------------------
summ = json.load(open(os.path.join(OUT, "neou_internal_cluster_summary.json")))
rows = list(csv.DictReader(open(os.path.join(OUT, "neou_internal_clusters.csv"))))
clab = {0: "c0 (SF1-lineage)", 1: "c3 (TPIT)", 2: "c4 (PIT1)"}
ccmap = {0: "0 (SF1-lineage)", 1: "3 (TPIT)", 2: "4 (PIT1)"}
st2 = []
for k in ["0", "1", "2"]:
    c = summ["cluster_characterization"][k]
    dom_lin = max(c["lineages"], key=c["lineages"].get)
    dom_stg = max(c["stages"], key=c["stages"].get)
    cross = c["cross_cohort_cluster_counts"]
    cross_str = ", ".join(f"c{k2}={v}" for k2, v in sorted(cross.items()) if v > 0)
    st2.append({
        "internal_cluster": f"int{k}",
        "n": c["n"],
        "nfpt_n": c["nfpt"],
        "nfpt_pct": round(100 * c["nfpt_frac"], 1),
        "mean_moller_key_beta": round(c["moller_key_beta"], 4),
        "mean_lineage_score": round(c["lineage_score"], 4),
        "dominant_lineage": dom_lin,
        "dominant_stage": dom_stg,
        "cross_cohort_centroid": ccmap[int(k)],
        "cross_cohort_counts": cross_str,
    })
with open(os.path.join(OUT, "supplementary_table2_internal_clusters.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(st2[0].keys()))
    w.writeheader(); w.writerows(st2)
print("[ST2] clusters:", [(r["internal_cluster"], r["n"], r["nfpt_n"], r["dominant_lineage"]) for r in st2])

# ---------------------------------------------------------------------------
# ST3 : Normalization sensitivity (index coding 1=c0,2=c2,3=c3,4=c4)
# ---------------------------------------------------------------------------
ns = list(csv.DictReader(open(os.path.join(G, "norm_sens_assignments.csv"))))
idx2cl = {"1": "c0", "2": "c2", "3": "c3", "4": "c4"}
methods = ["raw", "noob", "fun", "quant"]
mlab = {"raw": "preprocessRaw", "noob": "noob", "fun": "funnorm", "quant": "quantile"}
st3 = []
for m in methods:
    cnt = {f"c{i}": 0 for i in range(5)}
    for r in ns:
        cnt[idx2cl[r[m]]] += 1
    c0 = cnt["c0"]; n = len(ns)
    agree = round(100 * sum(1 for r in ns if r[m] == r["raw"]) / n, 1) if m != "raw" else 100.0
    st3.append({
        "normalization": mlab[m],
        "n_samples": n,
        "c0": cnt["c0"], "c2": cnt["c2"], "c3": cnt["c3"], "c4": cnt["c4"],
        "c0_fraction": round(c0 / n, 4),
        "agreement_vs_preprocessRaw_pct": agree,
    })
with open(os.path.join(OUT, "supplementary_table3_normalization.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(st3[0].keys()))
    w.writeheader(); w.writerows(st3)
# reassigned samples (flip under any norm vs preprocessRaw)
flips = [r["sample"] for r in ns if len({r["raw"], r["noob"], r["fun"], r["quant"]}) > 1]
print("[ST3] methods:", [(r["normalization"], r["c0"], r["c0_fraction"], r["agreement_vs_preprocessRaw_pct"]) for r in st3])
print("[ST3] reassigned samples (vs preprocessRaw):", flips)

# dump a combined json for embedding
summary = {
    "st1_pc1_ratio": round(float(ratio[0]), 4), "st1_pc2_ratio": round(float(ratio[1]), 4),
    "st1_cum20": round(float(ratio[:20].sum()), 4),
    "st2": st2, "st3": st3, "st3_reassigned": flips,
    "st2_best_k": summ["best_k"], "st2_silhouette_by_k": summ["silhouette_by_k"],
}
json.dump(summary, open(os.path.join(OUT, "supplementary_tables_summary.json"), "w"), indent=2)
print("[done] wrote 3 CSVs + summary json")
