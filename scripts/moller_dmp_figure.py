"""Figure: Møller 2025 562-DMP probe-level cross-cohort validation.

Panel A: direction-consistency heatmap of each Møller cluster (k1-k5) across the
         candidate groups of the three cohorts. Cell color = fraction of that
         cluster's DMPs whose beta-direction matches Møller's labeled trend.
         Significance (binomial p<0.05) marked with *.
Panel B: overall 562-DMP direction consistency per cohort with binomial p-values.

Outputs: results/Moller2025/moller_dmp_validation_figure.png (+ pdf)
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASE = "results/Moller2025"
mapping = json.load(open(f"{BASE}/moller_k_cluster_mapping.json"))
summary = json.load(open(f"{BASE}/moller_dmp_validation_summary.json"))

# --- Panel A matrix -------------------------------------------------------
rows = ["k1", "k2", "k3", "k4", "k5"]
# column = (display label, cohort key, group key)
cols = [
    ("GSE283928 c0", "GSE283928", "0"),
    ("GSE283928 c3", "GSE283928", "3"),
    ("GSE283928 c4", "GSE283928", "4"),
    ("GSE207937 NS", "GSE207937", "NS"),
    ("GSE207937 secr", "GSE207937", "secreting"),
    ("Neou int0", "Neou", "0"),
    ("Neou int1", "Neou", "1"),
    ("Neou int2", "Neou", "2"),
]

mat = np.full((len(rows), len(cols)), np.nan)
sig = np.zeros((len(rows), len(cols)), dtype=bool)
for i, k in enumerate(rows):
    for j, (_, coh, grp) in enumerate(cols):
        pg = mapping[coh][k]["per_group"].get(grp)
        if pg is None:
            continue
        mat[i, j] = pg["frac"]
        try:
            sig[i, j] = float(str(pg["p"])) < 0.05
        except Exception:
            sig[i, j] = False

# --- Build figure ---------------------------------------------------------
fig = plt.figure(figsize=(11.5, 5.6))
gs = fig.add_gridspec(1, 2, width_ratios=[2.05, 1.0], wspace=0.32)

# Panel A: heatmap
axA = fig.add_subplot(gs[0])
cmap = plt.cm.RdYlGn
im = axA.imshow(mat, cmap=cmap, vmin=0.0, vmax=1.0, aspect="auto")
axA.set_xticks(range(len(cols)))
axA.set_xticklabels([c[0] for c in cols], rotation=35, ha="right", fontsize=9)
axA.set_yticks(range(len(rows)))
axA.set_yticklabels(rows, fontsize=10)
axA.set_xlabel("Target cohort group", fontsize=10)
axA.set_ylabel("Møller 2025 cluster", fontsize=10)
axA.set_title("A. Probe-level direction consistency of 562 DMPs\n(fraction of DMPs matching Møller's labeled trend)", fontsize=10.5)

for i in range(len(rows)):
    for j in range(len(cols)):
        v = mat[i, j]
        if np.isnan(v):
            axA.add_patch(Rectangle((j-0.5, i-0.5), 1, 1, fill=True, color="0.85"))
            axA.text(j, i, "n/a", ha="center", va="center", fontsize=7, color="0.4")
            continue
        txt = f"{v:.2f}"
        if sig[i, j]:
            txt += "*"
        # choose text color by background luminance
        axA.text(j, i, txt, ha="center", va="center", fontsize=8.5,
                 color="black" if v < 0.55 or v > 0.85 else "black",
                 fontweight="bold" if sig[i, j] else "normal")

cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.02)
cb.set_label("direction consistency", fontsize=8.5)
axA.text(0.5, -0.27, "SF1-lineage Møller clusters (k1/k2/k3/k5) align with NFPT-enriched groups (c0/NS/int0);\n"
                     "k4 (HR 3.63, TPIT/PIT1, highest risk) aligns with secretory/TPIT groups (c3/c4/secr/int1). * p<0.05",
         transform=axA.transAxes, ha="center", fontsize=7.6, color="0.25")

# Panel B: overall per-cohort bars
axB = fig.add_subplot(gs[1])
cohorts = ["GSE283928", "GSE207937", "Neou"]
fracs = [summary["cohorts"][c]["frac_correct_direction"] for c in cohorts]
ps = [summary["cohorts"][c]["binom_p_all"] for c in cohorts]
labels = ["GSE283928\n(n=116)", "GSE207937\n(n=77)", "Neou\n(n=86)"]
colors = ["#2e7d32" if p < 0.05 else "#9e9e9e" for p in ps]
bars = axB.bar(range(3), fracs, color=colors, edgecolor="black", linewidth=0.8, width=0.6)
axB.axhline(0.5, color="crimson", linestyle="--", linewidth=1.0)
axB.text(2.45, 0.505, "chance = 0.5", color="crimson", fontsize=7.5, ha="right", va="bottom")
axB.set_xticks(range(3))
axB.set_xticklabels(labels, fontsize=8.5)
axB.set_ylim(0.45, 0.70)
axB.set_ylabel("Overall 562-DMP direction consistency", fontsize=9.5)
axB.set_title("B. Cross-cohort replication", fontsize=10.5)
for b, f, p in zip(bars, fracs, ps):
    pe = f"<1e-7" if p < 1e-7 else f"{p:.1e}"
    axB.text(b.get_x()+b.get_width()/2, f+0.008, f"{f:.3f}\n(p={pe})",
             ha="center", va="bottom", fontsize=7.6)

fig.suptitle("Probe-level validation of the Møller-2025 methylation taxonomy (562 DMPs) across three independent EPIC cohorts",
             fontsize=11.5, fontweight="bold", y=0.99)
fig.savefig(f"{BASE}/moller_dmp_validation_figure.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{BASE}/moller_dmp_validation_figure.pdf", bbox_inches="tight")
print("Saved figure:", f"{BASE}/moller_dmp_validation_figure.png")
print("Matrix:\n", np.round(mat, 3))
print("Significant cells:\n", sig)
