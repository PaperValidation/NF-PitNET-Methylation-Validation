import gzip, numpy as np
MAT = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/GSE283928_MatrixProcessed_Normalized.txt.gz"
with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
    hdr = f.readline().rstrip("\n").split("\t")
beta = [i for i, h in enumerate(hdr) if h.endswith("_Beta")]
dpv = [i for i, h in enumerate(hdr) if "Detection Pval" in h]
print("ncol total:", len(hdr) - 1, "| _Beta cols:", len(beta), "| _Detection Pval cols:", len(dpv))
print("first Beta:", hdr[beta[0]], "... last Beta:", hdr[beta[-1]])
print("first DPval:", hdr[dpv[0]], "... last DPval:", hdr[dpv[-1]])
rows = []; probes = []
with gzip.open(MAT, "rt", encoding="utf-8", errors="ignore") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        probes.append(p[0]); rows.append(p[1:])
X = np.array(rows, dtype=np.float64)
print("matrix shape:", X.shape)
B = X[:, [i - 1 for i in beta]]
D = X[:, [i - 1 for i in dpv]]
print("BETA  cols median", np.median(B), "min", B.min(), "max", B.max())
print("DPVAL cols median", np.median(D), "min", D.min(), "max", D.max())
print("DPVAL fraction <0.05:", np.mean(D < 0.05), "| fraction==0:", np.mean(D == 0))
for pr in ["cg00000109", "cg00000155"]:
    r = probes.index(pr)
    print(pr, "Beta[:4]=", np.round(B[r, :4], 3), "DPval[:4]=", np.round(D[r, :4], 4))
