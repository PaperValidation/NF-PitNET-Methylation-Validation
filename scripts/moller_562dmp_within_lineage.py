"""
Direction 1: recover Møller's within-lineage k-subdivision (k1-k5) from our
INDEPENDENTLY minfi-processed Møller beta (mendeley_beta_cg) using the 562
published cluster-specific DMPs (moller_dmp_probe_validation.csv).

Goal: the coarse c0-centroid test (§3.8) could not resolve the within-lineage
k-subdivision that lives inside the SF1-lineage c0 cluster. Here we test whether
the published 562-DMP coordinate system, applied to independently processed beta,
recovers Møller's original k-assignment. This is a cross-processing reproducibility
check of the taxonomy (same 117 samples, but a beta matrix derived entirely
independently of Møller's own normalization).

Two estimators, both leave-one-out (LOO) on the 117 samples:
  (A) nearest-centroid: per-k mean 562-DMP profile -> assign each left-out sample
      to the closest centroid (Euclidean). Transparent, hyperparameter-free.
  (B) random forest (sklearn) if available: 562 features, 5 classes, LOO predict.

Outputs:
  results/Moller2025/within_lineage_562dmp/moller_562dmp_assignment.csv
  results/Moller2025/within_lineage_562dmp/moller_562dmp_summary.json
  results/Moller2025/within_lineage_562dmp/moller_562dmp_confusion.png
"""
import os, json, glob
import numpy as np
import pandas as pd
from _paths import DATASET_PATH

BASE = r"C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
RES = os.path.join(BASE, "results/Moller2025/within_lineage_562dmp")
os.makedirs(RES, exist_ok=True)

BETA_GZ = os.path.join(DATASET_PATH, "Mendeley_Moller2025", "mendeley_beta_cg.csv.gz")
DMP_CSV = os.path.join(BASE, "results/Moller2025/moller_dmp_probe_validation.csv")
ASSIGN_CSV = os.path.join(BASE, "results/Moller2025/strict_projection/moller_strict_assignment.csv")

# ---------------------------------------------------------------- load 562 DMPs
dmp = pd.read_csv(DMP_CSV)
dmp_probes = dmp["probe"].astype(str).tolist()
n_dmp = len(dmp_probes)
assert n_dmp == 562, f"expected 562 DMPs, got {n_dmp}"
print(f"[DMP] loaded {n_dmp} cluster-specific DMPs")

# per-k membership + published direction (predicted_sign: +1 Hyper, -1 Hypo)
k_members = dmp.groupby("moller_cluster")["probe"].apply(list).to_dict()
print(f"[DMP] per-k membership: " + ", ".join(f"{k}:{len(v)}" for k,v in sorted(k_members.items())))

# ---------------------------------------------------------------- align k-labels
# The strict-projection assignment table already bridges beta IDAT basenames to
# Møller's Sample_ID and Kcluster_5. 6/117 samples lack a k-label + follow-up in
# Møller's supplement (NaN) and are dropped, matching the prognostic analysis
# (111 labeled samples).
sa = pd.read_csv(ASSIGN_CSV).dropna(subset=["Kcluster_5"])
VALID_K = {"k1", "k2", "k3", "k4", "k5"}
idat_to_k = {str(r["IDAT"]): str(r["Kcluster_5"])
             for _, r in sa.iterrows() if r["Kcluster_5"] in VALID_K}
print(f"[align] {len(idat_to_k)} IDAT->valid-k-label pairs (strict_assignment bridge)")

# ---------------------------------------------------------------- subset beta to 562 rows (chunked)
print(f"[beta] chunked read of {BETA_GZ} ...")
probe_set = set(dmp_probes)
chunks = []
total_rows = 0
reader = pd.read_csv(BETA_GZ, sep="\t", index_col=0, chunksize=40000)
for ch in reader:
    total_rows += len(ch)
    hit = ch.loc[ch.index.intersection(probe_set)]
    if len(hit):
        chunks.append(hit)
beta_sub = pd.concat(chunks)
beta_sub = beta_sub.loc[dmp_probes]   # enforce exact DMP order
print(f"[beta] scanned {total_rows} probes; kept {beta_sub.shape[0]} DMP rows x {beta_sub.shape[1]} samples")

# order samples by k-label availability (keep only those with a valid k-label)
cols = [c for c in beta_sub.columns if idat_to_k.get(c) in VALID_K]
print(f"[align] using {len(cols)} labeled samples (dropped {beta_sub.shape[1]-len(cols)} unlabeled)")
beta_sub = beta_sub[cols]
M = beta_sub.to_numpy(dtype=np.float32)   # 562 x n_labeled
sample_idat = np.array(cols)
sample_k = np.array([idat_to_k[c] for c in cols])
print(f"[beta] final matrix {M.shape}, k-distribution: " + ", ".join(f"{k}:{(sample_k==k).sum()}" for k in sorted(set(sample_k))))

# handle any NaN in DMP rows (impute column mean)
if np.isnan(M).any():
    col_ok = ~np.isnan(M)
    cm = np.nanmean(M, axis=1, keepdims=True)
    M = np.where(np.isnan(M), cm, M)

# ---------------------------------------------------------------- LOO nearest-centroid
K = sorted(set(sample_k))
kidx = {k: np.where(sample_k == k)[0] for k in K}
# centroids from full data (for descriptive), but LOO predictions exclude the test sample
pred_nc = np.empty(len(sample_k), dtype=object)
for i in range(len(sample_k)):
    tr = np.delete(np.arange(len(sample_k)), i)
    cents = np.stack([M[:, kidx[k]].mean(axis=1) for k in K])   # len(K) x 562
    d = np.linalg.norm(cents - M[:, i][None, :], axis=1)
    pred_nc[i] = K[int(np.argmin(d))]
acc_nc = float((pred_nc == sample_k).mean())
print(f"[NC-LOO] accuracy = {acc_nc:.3f}")

# ---------------------------------------------------------------- LOO random forest (if sklearn)
pred_rf = None; acc_rf = None; import_ok = False
try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import LeaveOneOut
    from sklearn.metrics import cohen_kappa_score, confusion_matrix, accuracy_score
    import_ok = True
    X = M.T   # 117 x 562
    y = sample_k
    pred_rf = np.empty(len(y), dtype=object)
    loo = LeaveOneOut()
    for tr, te in loo.split(X):
        clf = RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=-1)
        clf.fit(X[tr], y[tr])
        pred_rf[te] = clf.predict(X[te])[0]
    acc_rf = float(accuracy_score(y, pred_rf))
    kappa_rf = float(cohen_kappa_score(y, pred_rf))
    print(f"[RF-LOO] accuracy = {acc_rf:.3f}, kappa = {kappa_rf:.3f}")
except Exception as e:
    print(f"[RF] skipped ({e})")

# ---------------------------------------------------------------- confusion + agreement
def confusion(true, pred):
    cm = np.zeros((len(K), len(K)), dtype=int)
    for t, p in zip(true, pred):
        cm[K.index(t), K.index(p)] += 1
    return cm

cm_nc = confusion(sample_k, pred_nc)
# naive baseline: assign all to majority class
maj = max(set(sample_k), key=lambda k: (sample_k == k).sum())
baseline_acc = float((sample_k == maj).mean())

summary = {
    "n_samples": int(len(sample_k)),
    "n_dmp": int(n_dmp),
    "k_distribution": {k: int((sample_k == k).sum()) for k in K},
    "majority_class": maj,
    "naive_majority_accuracy": baseline_acc,
    "nearest_centroid_LOO_accuracy": acc_nc,
    "nearest_centroid_LOO_confusion": cm_nc.tolist(),
    "class_order": K,
}
if import_ok:
    cm_rf = confusion(sample_k, pred_rf)
    summary["random_forest_LOO_accuracy"] = acc_rf
    summary["random_forest_LOO_kappa"] = kappa_rf
    summary["random_forest_LOO_confusion"] = cm_rf.tolist()

# ---------------------------------------------------------------- write assignment csv
out = pd.DataFrame({
    "IDAT": sample_idat,
    "Sample_ID": [idat_to_k[c] for c in sample_idat],
    "actual_k": sample_k,
    "pred_k_centroid": pred_nc,
})
if import_ok:
    out["pred_k_rf"] = pred_rf
out.to_csv(os.path.join(RES, "moller_562dmp_assignment.csv"), index=False)

with open(os.path.join(RES, "moller_562dmp_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
print("[write] assignment + summary JSON")

# ---------------------------------------------------------------- confusion figure
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2 if import_ok else 1, figsize=(5*(2 if import_ok else 1), 4.2))
    def heat(ax, cm, title):
        im = ax.imshow(cm, cmap="Blues")
        ax.set_xticks(range(len(K))); ax.set_xticklabels(K)
        ax.set_yticks(range(len(K))); ax.set_yticklabels(K)
        ax.set_xlabel("Predicted k"); ax.set_ylabel("Actual Møller k")
        ax.set_title(title)
        for r in range(len(K)):
            for c in range(len(K)):
                ax.text(c, r, cm[r, c], ha="center", va="center",
                        color="white" if cm[r, c] > cm.max()/2 else "black", fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    heat(ax[0] if import_ok else ax, cm_nc, f"Nearest-centroid LOO (acc={acc_nc:.0%})")
    if import_ok:
        heat(ax[1], cm_rf, f"RandomForest LOO (acc={acc_rf:.0%})")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "moller_562dmp_confusion.png"), dpi=150)
    print("[write] confusion figure")
except Exception as e:
    print(f"[fig] skipped ({e})")

print("[done]")
