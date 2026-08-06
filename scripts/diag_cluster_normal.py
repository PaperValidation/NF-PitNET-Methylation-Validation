import gzip, json, os
import numpy as np
import pandas as pd

BASE = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
MATP = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
CL   = os.path.join(BASE, "results/GSE283928/cluster_assignment_k5.csv")
OUT  = os.path.join(BASE, "results/GSE283928/cluster_diag.json")

cl = pd.read_csv(CL)
cl_map = dict(zip(cl["sample_id"], cl["cluster"]))
clusters = sorted(set(cl_map.values()))
print("clusters:", {c: sum(1 for v in cl_map.values() if v == c) for c in clusters})

samples = None
# We compute, per cluster, global median beta and % of probes with beta<0.2 (open)
# streaming: accumulate per-cluster order statistics is hard; instead keep full matrix small enough? 727k*236*4 = 1.37GB
# Use subsample of first 100k probes to estimate global median per cluster reliably (median is robust to sampling)
cluster_below = {c: 0.0 for c in clusters}
cluster_count = {c: 0.0 for c in clusters}
cluster_vals_sum = {c: 0.0 for c in clusters}
import random
rng = random.Random(42)
SUB = 120000
kept = 0
with gzip.open(MATP, "rt", encoding="utf-8", errors="ignore") as f:
    header = f.readline().rstrip("\n").split("\t")
    samples = header[1:]
    arr = None
    for i, line in enumerate(f):
        if kept >= SUB:
            break
        parts = line.rstrip("\n").split("\t")
        vals = np.array(parts[1:], dtype=np.float32)
        # assign to clusters
        for c in clusters:
            msk = np.array([cl_map[s] == c for s in samples])
            cluster_vals_sum[c] += vals[msk].sum()
            cluster_count[c] += msk.sum()
            cluster_below[c] += (vals[msk] < 0.2).sum()
        kept += 1

res = {}
for c in clusters:
    n_probe_sample = cluster_count[c]
    res[int(c)] = {
        "n_samples": int(sum(1 for v in cl_map.values() if v == c)),
        "global_mean_beta_sub": round(float(cluster_vals_sum[c] / n_probe_sample), 4) if n_probe_sample else None,
        "frac_beta_lt_0.2_sub": round(float(cluster_below[c] / n_probe_sample), 4) if n_probe_sample else None
    }
with open(OUT, "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
