import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os, json

from _paths import DATASET_PATH as BASE
OUT = os.path.join(BASE, "results", "overview_figure.png")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

fig = plt.figure(figsize=(14, 9))
fig.suptitle("NF-PitNET Post-surgical Regrowth — Multi-omics Surrogate Framework\n(3 major milestones completed)", fontsize=15, fontweight="bold")

# Panel A: GSE283928 cluster sizes + normal identification
axA = fig.add_subplot(2, 3, 1)
clusters = ["0", "1 (normal)", "2", "3", "4"]
sizes = [70, 118, 2, 22, 24]
colors = ["#4C72B0", "#C44E52", "#55A868", "#8172B3", "#CCB974"]
bars = axA.bar(clusters, sizes, color=colors)
axA.set_title("A. GSE283928 k=5 clusters (Belakhoua 2025)\ncluster 1 = normal control (beta~0)", fontsize=11)
axA.set_ylabel("samples")
for b, s in zip(bars, sizes):
    axA.text(b.get_x()+b.get_width()/2, s+1, str(s), ha="center", fontsize=9)

# Panel B: lineage dominance score by tumor cluster
axB = fig.add_subplot(2, 3, 2)
ld = json.load(open(os.path.join(BASE, "results/GSE283928/lineage_cluster_validation.json")))
sd = ld["SF1_dominance_score_by_cluster"]
tumors = [k for k in sd if k != "1"]
vals = [sd[k]["mean_SF1_score"] for k in tumors]
axB.bar([f"c{k}" for k in tumors], vals, color="#4C72B0")
axB.axhline(0, color="k", lw=0.8)
axB.set_title("B. SF1-dominance score by tumor cluster\n(mixed lineage = Møller 4 SF1 + 1 non-SF1)", fontsize=11)
axB.set_ylabel("SF1 - (TPIT+PIT1)/2 methylation score")

# Panel C: Møller key genes overlap with GSE209903
axC = fig.add_subplot(2, 3, 3)
moller = 13
sig_in = 10
axC.bar(["Møller key genes", "Significant in GSE209903"], [moller, sig_in], color=["#C44E52", "#55A868"])
axC.set_title("C. Møller (2025) signature anchors\n10/13 key genes significant in GSE209903", fontsize=11)
axC.set_ylabel("gene count")
for i, v in enumerate([moller, sig_in]):
    axC.text(i, v+0.2, str(v), ha="center", fontsize=10)

# Panel D: drug repurposing priority
axD = fig.add_subplot(2, 3, 4)
drugs = ["Palbociclib\n(CDK4/6)", "Trametinib\n(MEK)", "Anti-SPP1/\nITGAV-ITGB1", "Linifanib\n(VEGFR)", "Eplerenone\n(MR)"]
prio = [3, 3, 3, 2, 2]
barcol = ["#C44E52" if p==3 else "#CCB974" for p in prio]
axD.barh(drugs[::-1], prio[::-1], color=barcol[::-1])
axD.set_xlim(0, 3.4)
axD.set_title("D. Drug repurposing candidates\n(★★★ highest evidence)", fontsize=11)
axD.set_xlabel("priority")

# Panel E: cross-omics methylation of NFPT-up program
axE = fig.add_subplot(2, 3, 5)
cm = json.load(open(os.path.join(BASE, "results/risk_model/nfpt_program_by_cluster.json")))
ents = ["TUMOR_all", "NORMAL_all"]
mv = [cm[e]["mean_methylation_NFPTup_program"] for e in ents]
axE.bar(["Tumor (n=118)", "Normal (n=118)"], mv, color=["#4C72B0", "#C44E52"])
axE.set_title("E. NFPT-up program methylation\n(GSE209903 genes -> EPIC probes)", fontsize=11)
axE.set_ylabel("mean beta")
for i, v in enumerate(mv):
    axE.text(i, v+0.01, f"{v:.2f}", ha="center", fontsize=9)

# Panel F: data cohorts summary
axF = fig.add_subplot(2, 3, 6)
cohorts = ["GSE51618\n(invasion)", "GSE209903\n(NFPT RNA)", "GSE207937\n(methyl)", "GSE283928\n(methyl)"]
n = [10, 65, 77, 236]
axF.bar(cohorts, n, color="#8172B3")
axF.set_title("F. Open cohorts acquired\n(874MB + 490MB methyl matrices)", fontsize=11)
axF.set_ylabel("samples")
for i, v in enumerate(n):
    axF.text(i, v+3, str(v), ha="center", fontsize=9)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(OUT, dpi=130, bbox_inches="tight")
print("Saved", OUT)
