#!/usr/bin/env python3
"""External lineage validation of the c0 (SF1-lineage) mapping on an INDEPENDENT
NFPA-only methylation cohort: GSE115783 (Bujko/Kober, Illumina HumanMethylation450,
34 non-functioning pituitary adenomas + 6 healthy pituitary controls).

We project every GSE115783 sample onto the GSE283928-derived c0/c2/c3/c4 centroids
(nearest-centroid, NO k-labels) and ask: do independent NFPA fall predominantly in
c0, replicating the 86% c0 finding from the Møller 2025 raw-idat strict test?

The 6 healthy-pituitary samples serve as a built-in specificity control: they should
NOT map to c0 (a tumor SF1-lineage identity).

Key facts:
- GSE283928 centroids are EPIC 850K; GSE115783 is HM450. Intersection ~450k cg.
- GEO series matrix is TAB-separated (R/write.table export); read with sep="\t".
- Beta values are normalized [0,1] on both sides (submitter ChAMP-processed).
"""
import os, json, gzip, re
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from _paths import DATASET_PATH

BASE = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
CENT = os.path.join(DATASET_PATH, "Moller2025_derived", "gse283928_centroids.npz")
GSE  = os.path.join(DATASET_PATH, "GSE115783", "GSE115783_series_matrix.txt.gz")
OUTD = os.path.join(BASE, "results/Moller2025/gse115783_validation")
os.makedirs(OUTD, exist_ok=True)

CLUSTERS = ["c0", "c2", "c3", "c4"]   # c1 is the detection-pvalue artifact, excluded

# ---------------------------------------------------------------- load centroids
z = np.load(CENT, allow_pickle=True)
cent_probes = z["probes"].astype(str)
C = {c: z[f"centroid_{c}"].astype(np.float32) for c in CLUSTERS}
cent_idx = {p: i for i, p in enumerate(cent_probes)}
print(f"[centroid] probes={len(cent_probes)}  clusters={CLUSTERS}")

# ---------------------------------------------------------------- parse annotations (tissue)
tissue_line = None
table_begin = None
with gzip.open(GSE, "rt", encoding="utf-8", errors="replace") as f:
    for i, line in enumerate(f):
        if line.startswith("!Sample_characteristics_ch1"):
            fld = line.rstrip("\n").split("\t")
            vals = [v.strip().strip('"') for v in fld[1:]]
            if vals and vals[0].lower().startswith("tissue"):
                tissue_line = vals
        if line.startswith("!series_matrix_table_begin"):
            table_begin = i
            break
assert table_begin is not None, "series_matrix_table_begin marker not found"
print(f"[annot] table_begin at line {table_begin}; tissue labels captured: {len(tissue_line) if tissue_line else 0}")

# ---------------------------------------------------------------- load beta (TSV!)
beta = pd.read_csv(GSE, sep="\t", skiprows=table_begin + 1, index_col=0, on_bad_lines="skip")
beta = beta.astype(np.float32)
print(f"[gse] raw matrix {beta.shape}  (rows=probes, cols=samples)")

# ---------------------------------------------------------------- classify columns
def classify(label):
    if label is None:
        return "unknown"
    l = label.lower()
    if "healthy" in l or "normal" in l:
        return "normal"
    if "nonfunctioning" in l or "adenoma" in l or "nfpa" in l:
        return "tumor"
    return "other"

groups = {}
if tissue_line is not None:
    assert len(tissue_line) == beta.shape[1], \
        f"tissue labels ({len(tissue_line)}) != columns ({beta.shape[1]})"
    for col, lab in zip(beta.columns, tissue_line):
        groups[col] = classify(lab)
else:
    for col in beta.columns:
        groups[col] = "unknown"
print("[annot] group counts:", dict(Counter(groups.values())))

v = beta.to_numpy()
print(f"[gse] beta range [{np.nanmin(v):.3f}, {np.nanmax(v):.3f}]")

# ---------------------------------------------------------------- common probes
beta_probe_set = set(beta.index)
common = [p for p in cent_probes if p in beta_probe_set]
idx_c = np.array([cent_idx[p] for p in common])               # centroid row pos
idx_g = np.array([beta.index.get_loc(p) for p in common])     # beta row pos
print(f"[common] {len(common)} shared cg probes (centroid EPIC ∩ GSE115783 HM450)")

Cc = np.stack([C[c][idx_c] for c in CLUSTERS], axis=1)        # (P, 4) centroids
Mg = beta.iloc[idx_g, :].to_numpy(np.float32)                 # (P, N)

# ---------------------------------------------------------------- NaN hygiene
row_ok = ~np.isnan(Mg).all(axis=1)
Mg = Mg[row_ok]
Cc = Cc[row_ok]
col_mean = np.nanmean(Mg, axis=1, keepdims=True)
Mg = np.where(np.isnan(Mg), col_mean, Mg)
print(f"[clean] after NaN drop/impute -> {Mg.shape}")

# ---------------------------------------------------------------- nearest centroid
d2 = ((Mg[:, :, None] - Cc[:, None, :]) ** 2).sum(axis=0)     # (N, 4)
assign = d2.argmin(axis=1)
dist = d2.min(axis=1)
sample_names = np.array(beta.columns.tolist())
assign_labels = np.array(CLUSTERS)[assign]

# ---------------------------------------------------------------- stratified report
assign_df = pd.DataFrame({"Sample": sample_names,
                          "group": [groups[s] for s in sample_names],
                          "assigned_cluster": assign_labels,
                          "min_distance": dist})

summary = {"cohort": "GSE115783",
           "platform": "Illumina HumanMethylation450 (HM450)",
           "tumor_type": "non-functioning pituitary adenoma (NFPA) + healthy pituitary control",
           "centroid_source": "GSE283928 (Belakhoua 2025, EPIC 850K)",
           "common_probes": int(len(common)),
           "reference_moller_c0_fraction": 0.86,
           "n_total": int(len(sample_names)),
           "groups": {}}
for g in ["tumor", "normal", "other", "unknown"]:
    sub = assign_df[assign_df.group == g]
    if len(sub) == 0:
        continue
    counts = {c: int((sub.assigned_cluster == c).sum()) for c in CLUSTERS}
    summary["groups"][g] = {"n": int(len(sub)),
                            "assignment_counts": counts,
                            "c0_fraction": round(counts["c0"] / len(sub), 4)}
counts_all = {c: int((assign_labels == c).sum()) for c in CLUSTERS}
summary["overall"] = {"n": int(len(assign_labels)),
                      "assignment_counts": counts_all,
                      "c0_fraction": round(counts_all["c0"] / len(assign_labels), 4)}

print("\n=== GSE115783 lineage validation ===")
print(json.dumps(summary, indent=2))

# ---------------------------------------------------------------- outputs
assign_df.to_csv(os.path.join(OUTD, "gse115783_assignment.csv"), index=False)
with open(os.path.join(OUTD, "gse115783_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)

# bar chart: tumor vs normal assignment profile
fig, ax = plt.subplots(figsize=(7, 4.5))
xs = np.arange(len(CLUSTERS))
palette = {"tumor": "#2c7fb8", "normal": "#d95f0e", "other": "#999999", "unknown": "#cccccc"}
gi = 0
offsets = {"tumor": -0.19, "normal": 0.19}
for g in ["tumor", "normal", "other", "unknown"]:
    if g not in summary["groups"]:
        continue
    cnt = summary["groups"][g]["assignment_counts"]
    ax.bar(xs + offsets.get(g, 0), [cnt[c] for c in CLUSTERS], width=0.38,
           label=f"{g} (n={summary['groups'][g]['n']})", color=palette[g])
    gi += 1
ax.set_xticks(xs)
ax.set_xticklabels(CLUSTERS)
ax.set_xlabel("Nearest centroid")
ax.set_ylabel("Number of samples")
tc0 = summary["groups"].get("tumor", {}).get("c0_fraction", 0) * 100
nc0 = summary["groups"].get("normal", {}).get("c0_fraction", 0) * 100
ax.set_title(f"GSE115783 NFPA -> GSE283928 centroids\nTumor c0 = {tc0:.0f}%  |  Normal c0 = {nc0:.0f}%")
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTD, "gse115783_lineage_bar.png"), dpi=150)
plt.close()
print(f"[done] outputs in {OUTD}")
