import gzip, json, os
import numpy as np
import pandas as pd

BASE = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
MATP = os.path.join(BASE, "GSE283928_MatrixProcessed_Normalized.txt.gz")
PROBES = os.path.join(BASE, "results/risk_model/nfpt_up_epic_probes.csv")
CL   = os.path.join(BASE, "results/GSE283928/cluster_assignment_k5.csv")
OUT  = os.path.join(BASE, "results/risk_model/nfpt_program_by_cluster.json")

probes_df = pd.read_csv(PROBES)
probe_list = probes_df.Probe.tolist()
pidx = {p: i for i, p in enumerate(probe_list)}

cl = pd.read_csv(CL)
cl_map = dict(zip(cl["sample_id"], cl["cluster"]))
# tumor clusters: exclude cluster 1 (normal control)
samples = None
mat = np.zeros((len(probe_list), 0), dtype=np.float32)
with gzip.open(MATP, "rt", encoding="utf-8", errors="ignore") as f:
    header = f.readline().rstrip("\n").split("\t")
    samples = header[1:]
    mat = np.zeros((len(probe_list), len(samples)), dtype=np.float32)
    for line in f:
        parts = line.rstrip("\n").split("\t")
        pid = parts[0]
        if pid in pidx:
            mat[pidx[pid], :] = np.array(parts[1:], dtype=np.float32)

cluster_of = np.array([cl_map[s] for s in samples])
# For each TUMOR cluster, mean methylation of NFPT-up program probes
res = {}
for c in sorted(set(cluster_of)):
    msk = cluster_of == c
    n_s = int(msk.sum())
    mean_beta = float(mat[:, msk].mean()) if n_s > 0 else None
    res[int(c)] = {"n_samples": n_s, "mean_methylation_NFPTup_program": round(mean_beta, 4) if mean_beta else None}

# Also tumor-only (exclude cluster1) overall vs normal (cluster1)
tumor_msk = cluster_of != 1
normal_msk = cluster_of == 1
res["TUMOR_all"] = {"n_samples": int(tumor_msk.sum()),
                    "mean_methylation_NFPTup_program": round(float(mat[:, tumor_msk].mean()), 4)}
res["NORMAL_all"] = {"n_samples": int(normal_msk.sum()),
                     "mean_methylation_NFPTup_program": round(float(mat[:, normal_msk].mean()), 4)}

with open(OUT, "w") as fh:
    json.dump(res, fh, indent=2)
print(json.dumps(res, indent=2))
print("\nInterpretation: NFPT-up transcriptional program -> if HYPOMETHYLATED (low beta) in tumor vs normal,")
print("it implies these loci are open/active in tumors, consistent with NFPT-up expression. Higher methylation = silenced.")
