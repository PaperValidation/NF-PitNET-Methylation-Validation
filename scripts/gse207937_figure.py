#!/usr/bin/env python3
"""Figure for GSE207937 NFPT methylome validation (4 panels)."""
import gzip, csv, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from _paths import DATASET_PATH as BASE
MAT = os.path.join(BASE, "GSE207937_series_matrix.txt.gz")
OUT = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/results/GSE207937"

def clean(s): return s.strip().strip('"')
def m2beta(m): return 1.0/(1.0+np.power(2.0,-m))

# parse matrix (beta)
header=None; probes=[]; rows=[]
with gzip.open(MAT,'rt',encoding='utf-8',errors='ignore') as f:
    for line in f:
        if line.startswith('!') or line.startswith('#') or line.strip()=='':
            continue
        parts=[clean(p) for p in line.rstrip('\n').split('\t')]
        if header is None:
            header=parts; continue
        probes.append(parts[0]); rows.append(parts[1:])
Xb=m2beta(np.array(rows,dtype=np.float32))
sample_ids=header[1:]

# labels from saved csv
subtypes={}; cluster={}; moller={}
with open(os.path.join(OUT,"sample_labels.csv")) as f:
    for r in csv.DictReader(f):
        subtypes[r["sample_id"]]=r["subtype"]
        cluster[r["sample_id"]]=int(r["cluster"])
        moller[r["sample_id"]]=float(r["moller_key_beta"])
sub = np.array([subtypes[s] for s in sample_ids])
cl  = np.array([cluster[s] for s in sample_ids])
mol = np.array([moller[s] for s in sample_ids])

# PCA top2
var=Xb.var(axis=1); k=int(len(var)*0.10); thr=np.partition(var,-k)[-k]
Xf=Xb[var>=thr,:]; Xc=Xf-Xf.mean(axis=0,keepdims=True)
pca=PCA(n_components=2,random_state=42); PC=pca.fit_transform(Xc.T)
evr=pca.explained_variance_ratio_

# color map
cmap={"NS":"#d1495b","GH":"#30638e","ACTH":"#edae49","NA":"#999999"}
cols=[cmap.get(s,"#999999") for s in sub]

fig,axes=plt.subplots(2,2,figsize=(11,9))
fig.suptitle("GSE207937 (Silva-Junior 2023) methylome stratification & NFPT validation",
             fontsize=13,fontweight='bold')

# Panel 1: PCA by subtype
ax=axes[0,0]
for st in ["NS","GH","ACTH","NA"]:
    m=sub==st
    ax.scatter(PC[m,0],PC[m,1],s=42,label=f"{st} (n={m.sum()})",
               color=cmap[st],edgecolor='k',linewidth=0.3,alpha=0.85)
ax.set_xlabel(f"PC1 ({evr[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({evr[1]*100:.1f}%)")
ax.set_title("A. PCA by subtype (NFPT=NS)"); ax.legend(fontsize=8,frameon=False)
ax.grid(alpha=0.25)

# Panel 2: cluster x subtype stacked bar
ax=axes[0,1]
order=["GH","ACTH","NS","NA"]; cl_order=[0,1,2,3,4]
bottom=np.zeros(len(cl_order))
for st in order:
    vals=[int(((cl==c)&(sub==st)).sum()) for c in cl_order]
    ax.bar([str(c) for c in cl_order],vals,bottom=bottom,label=st,color=cmap[st])
    bottom+=np.array(vals)
ax.set_xlabel("Unsupervised cluster (k=5)"); ax.set_ylabel("n samples")
ax.set_title("B. Cluster x subtype composition"); ax.legend(fontsize=8,frameon=False)

# Panel 3: Moller-key methyl by subtype
ax=axes[1,0]
mk={}
with open(os.path.join(OUT,"moller_key_by_subtype.json")) as f:
    mk=json.load(f)
sts=["GH","ACTH","NS"]; means=[mk[s] for s in sts]
bars=ax.bar(sts,means,color=[cmap[s] for s in sts],edgecolor='k',linewidth=0.5)
for b,v in zip(bars,means): ax.text(b.get_x()+b.get_width()/2,v+0.002,f"{v:.3f}",ha='center',fontsize=9)
ax.set_ylabel("Mean beta @ Moller-2025 key genes")
ax.set_title("C. Moller key-gene methylation by subtype")
ax.set_ylim(0.55,0.65)

# Panel 4: cross-cohort assignment heatmap
ax=axes[1,1]
xtab={}
with open(os.path.join(OUT,"crosscohort_assignment.csv")) as f:
    for r in csv.DictReader(f):
        xtab.setdefault(r["subtype"],{})
        # count assigned euclidean cluster
        c=r["assigned_cluster_eucl"]
        xtab[r["subtype"]][c]=xtab[r["subtype"]].get(c,0)+1
clusters=[0,2,3,4]
mat=np.array([[xtab.get(st,{}).get(str(c),0) for c in clusters] for st in ["GH","ACTH","NS"]],
             dtype=float)
# row-normalize for visualization of preference
rown=mat.sum(axis=1,keepdims=True); rown[rown==0]=1
matn=mat/rown
im=ax.imshow(matn,cmap='YlOrRd',aspect='auto',vmin=0,vmax=1)
ax.set_xticks(range(len(clusters))); ax.set_xticklabels([f"GSE283928 c{c}" for c in clusters])
ax.set_yticks(range(3)); ax.set_yticklabels(["GH","ACTH","NS(NFPT)"])
ax.set_title("D. Cross-cohort projection (row-normalized)")
for i in range(3):
    for j in range(len(clusters)):
        ax.text(j,i,f"{int(mat[i,j])}",ha='center',va='center',fontsize=9,
                color='black' if matn[i,j]<0.6 else 'white')
fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04,label="fraction within subtype")

plt.tight_layout(rect=[0,0,1,0.97])
fig.savefig(os.path.join(OUT,"gse207937_figure.png"),dpi=150,bbox_inches='tight')
print("saved", os.path.join(OUT,"gse207937_figure.png"))
