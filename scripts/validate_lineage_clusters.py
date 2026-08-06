import gzip, csv, json, os
import numpy as np
import pandas as pd

BASE = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
MATP = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
MK   = os.path.join(BASE, "results/GSE283928/lineage_marker_probes.csv")
CL   = os.path.join(BASE, "results/GSE283928/cluster_assignment_k5.csv")
OUT  = os.path.join(BASE, "results/GSE283928/lineage_cluster_validation.json")

print("Loading marker probe sets ...")
mk = pd.read_csv(MK)
sets = {s: mk[mk.set == s].probe.tolist() for s in mk.set.unique()}
for s, p in sets.items():
    print(f"  {s}: {len(p)} probes")

print("Loading cluster assignment ...")
cl = pd.read_csv(CL)  # columns: sample, cluster
cl_map = dict(zip(cl["sample_id"], cl["cluster"]))

# stream the big matrix, keep only marker probes (union)
all_probes = sorted(set(sum(sets.values(), [])))
probe_idx = {p: i for i, p in enumerate(all_probes)}
print(f"Total unique marker probes: {len(all_probes)}")

samples = None
betas = np.zeros((len(all_probes), 0), dtype=np.float32)
with gzip.open(MATP, "rt", encoding="utf-8", errors="ignore") as f:
    header = f.readline().rstrip("\n").split("\t")
    samples = header[1:]
    print(f"Matrix samples: {len(samples)}")
    # preallocate
    mat = np.zeros((len(all_probes), len(samples)), dtype=np.float32)
    for i, line in enumerate(f):
        parts = line.rstrip("\n").split("\t")
        pid = parts[0]
        if pid in probe_idx:
            mat[probe_idx[pid], :] = np.array(parts[1:], dtype=np.float32)
    print("Matrix slice built:", mat.shape)

# mean methylation per set per sample
cluster_of_sample = np.array([cl_map.get(s, -1) for s in samples])
result = {}
for s, plist in sets.items():
    idx = [probe_idx[p] for p in plist if p in probe_idx]
    mean_beta = mat[idx, :].mean(axis=0)  # per sample
    # per cluster mean
    clusters = sorted(set(cluster_of_sample))
    per_cluster = {}
    for c in clusters:
        msk = cluster_of_sample == c
        if msk.sum() > 0:
            per_cluster[int(c)] = round(float(mean_beta[msk].mean()), 4)
    result[s] = {
        "n_probes": len(idx),
        "per_cluster_mean_beta": per_cluster,
        "global_mean": round(float(mean_beta.mean()), 4)
    }

# Also compute a "lineage score" per sample: SF1 - (TPIT+PIT1)/2
sf1 = sets["SF1_lineage"]; tpit = sets["TPIT_lineage"]; pit1 = sets["PIT1_lineage"]
sf1_idx = [probe_idx[p] for p in sf1 if p in probe_idx]
tp_idx  = [probe_idx[p] for p in tpit if p in probe_idx]
pi_idx  = [probe_idx[p] for p in pit1 if p in probe_idx]
sf1_b = mat[sf1_idx,:].mean(0)
tp_b  = mat[tp_idx,:].mean(0)
pi_b  = mat[pi_idx,:].mean(0)
lineage_score = sf1_b - (tp_b + pi_b)/2.0  # high=SF1 dominant
# map to clusters
cluster_lineage = {}
for c in sorted(set(cluster_of_sample)):
    msk = cluster_of_sample == c
    if msk.sum() > 0:
        cluster_lineage[int(c)] = {
            "n": int(msk.sum()),
            "mean_SF1_score": round(float(lineage_score[msk].mean()), 4)
        }

result["SF1_dominance_score_by_cluster"] = cluster_lineage

with open(OUT, "w") as fh:
    json.dump(result, fh, indent=2)
print("Saved", OUT)
print(json.dumps(result, indent=2))
