import gzip, json, os
import numpy as np
import pandas as pd

from _paths import DATASET_PATH as BASE, DATASET_PATH_TEMP as TEMP
MATP = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
CL   = os.path.join(TEMP, "results/GSE283928/cluster_assignment_k5.csv")
OUT  = os.path.join(TEMP, "results/GSE283928/cluster_global_mean_beta.json")

cl = pd.read_csv(CL)
cl_map = dict(zip(cl["sample_id"], cl["cluster"]))

samples = None
# global mean beta per cluster (stream)
cluster_sum = {}
cluster_cnt = {}
with gzip.open(MATP, "rt", encoding="utf-8", errors="ignore") as f:
    header = f.readline().rstrip("\n").split("\t")
    samples = header[1:]
    for line in f:
        parts = line.rstrip("\n").split("\t")
        vals = np.array(parts[1:], dtype=np.float32)
        for c in set(cl_map[s] for s in samples):
            msk = np.array([cl_map[s] == c for s in samples])
            if c not in cluster_sum:
                cluster_sum[c] = np.zeros(len(samples), dtype=np.float64)
            cluster_sum[c] += vals
            cluster_cnt[c] = msk.sum()

res = {}
for c in sorted(cluster_sum):
    res[int(c)] = {
        "n": int(cluster_cnt[c]),
        "global_mean_beta": round(float(cluster_sum[c].sum() / (cluster_cnt[c] * len(samples))), 4)
    }
with open(OUT, "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
