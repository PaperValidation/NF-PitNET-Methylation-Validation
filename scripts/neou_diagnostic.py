import csv, numpy as np, os
NE="C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
OUT=os.path.join(NE,"results","Neou")
# distance profile
rows=list(csv.DictReader(open(os.path.join(OUT,"neou_crosscohort_assignment.csv"))))
clusters=[c for c in rows[0].keys() if c.startswith("dist_c")]
prof={c:[] for c in clusters}
prof_nfpt={c:[] for c in clusters}
prof_func={c:[] for c in clusters}
for r in rows:
    nfpt=int(r["nfpt"])
    for c in clusters:
        v=float(r[c]); prof[c].append(v)
        (prof_nfpt if nfpt==1 else prof_func)[c].append(v)
print("=== Neou mean distance to GSE283928 tumor centroids ===")
for c in clusters:
    a=np.array(prof[c]); an=np.array(prof_nfpt[c]); af=np.array(prof_func[c])
    print(f"  {c}: all={a.mean():.3f}+-{a.std():.3f} | NFPT={an.mean():.3f} | Func={af.mean():.3f}")
mn=min(np.mean(prof[c]) for c in clusters)
print(f"  -> min-mean cluster (all Neou): {[c for c in clusters if abs(np.mean(prof[c])-mn)<1e-9]}")

# PCA cohort positions
pc=list(csv.DictReader(open(os.path.join(OUT,"joint_pca_coords.csv"))))
from collections import defaultdict
co=defaultdict(list)
for r in pc:
    co[r["cohort"]].append([float(r[f"PC{i+1}"]) for i in range(20)])
print("\n=== Joint PCA cohort centroids (PC1,PC2) + spread ===")
for k,v in co.items():
    v=np.array(v)
    print(f"  {k}: n={v.shape[0]} meanPC1={v[:,0].mean():.3f} meanPC2={v[:,1].mean():.3f} | stdPC1={v[:,0].std():.3f} stdPC2={v[:,1].std():.3f}")
# overlap: distance between cohort centroids in PC1-PC2
keys=list(co.keys())
import itertools
print("  pairwise centroid dist (PC1-PC2 Euclidean):")
for a,b in itertools.combinations(keys,2):
    va=np.array(co[a]).mean(0)[:2]; vb=np.array(co[b]).mean(0)[:2]
    print(f"    {a}-{b}: {np.linalg.norm(va-vb):.3f}")
