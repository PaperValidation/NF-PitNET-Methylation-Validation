#!/usr/bin/env python3
"""Figure for Neou 2020 cross-cohort methylation validation (milestone 7).
4 panels:
  (A) Joint 3-cohort PCA scatter (PC1 vs PC2) - cohort-colored, NFPT shaped
  (B) Cross-cohort projection heatmap: per-sample similarity to GSE283928 4 tumor centroids
  (C) Moller-2025 key-gene methylation (mean +/- SD): Neou NFPT/Func vs GSE207937 NS/other
  (D) Lineage score boxplot (same 4 groups)
"""
import csv, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NE = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
OUT = os.path.join(NE, "results", "Neou")
G207_LAB = os.path.join(NE, "results", "GSE207937", "sample_labels.csv")

# ---------- load ----------
def rdcsv(p):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f))

neu_ass = rdcsv(os.path.join(OUT, "neou_crosscohort_assignment.csv"))
neu_lin = rdcsv(os.path.join(OUT, "neou_lineage_scores.csv"))
g207_ass = rdcsv(os.path.join(OUT, "gse207937_crosscohort_assignment.csv"))
g207_lab = rdcsv(G207_LAB)
pc = rdcsv(os.path.join(OUT, "joint_pca_coords.csv"))
clusters = ["0", "2", "3", "4"]

# Neou per-sample
neu_nfpt = np.array([int(r["nfpt"]) for r in neu_lin])
neu_mol = np.array([float(r["moller_key_beta"]) for r in neu_lin])
neu_ls = np.array([float(r["lineage_score"]) for r in neu_lin])
# GSE207937 per-sample
g207_sub = [r["subtype"] for r in g207_lab]
g207_nfpt = np.array([int(r["nfpt"]) for r in g207_lab])
g207_mol = np.array([float(r["moller_key_beta"]) for r in g207_lab])
g207_ls = np.array([float(r["lineage_score"]) for r in g207_lab])

# groups
groups = [
    ("Neou\nNFPT", neu_mol[neu_nfpt == 1], neu_ls[neu_nfpt == 1]),
    ("Neou\nFunctioning", neu_mol[neu_nfpt == 0], neu_ls[neu_nfpt == 0]),
    ("GSE207937\nNS (NFPT)", g207_mol[g207_nfpt == 1], g207_ls[g207_nfpt == 1]),
    ("GSE207937\nother", g207_mol[g207_nfpt == 0], g207_ls[g207_nfpt == 0]),
]

# ---------- colors ----------
C_NEU = "#d62728"; C_283 = "#1f77b4"; C_207 = "#2ca02c"

# ============ FIGURE ============
fig, ax = plt.subplots(2, 2, figsize=(13, 11))
fig.suptitle("Neou 2020 (E-MTAB-7762) cross-cohort methylation validation\n"
             "Third independent EPIC 850K cohort projected onto the MOLLER/BELAKHOUA taxonomy",
             fontsize=13, fontweight="bold")

# (A) PCA scatter
coh = np.array([r["cohort"] for r in pc])
nf = np.array([int(r["nfpt"]) for r in pc])
PC1 = np.array([float(r["PC1"]) for r in pc])
PC2 = np.array([float(r["PC2"]) for r in pc])
for c, col in [("Neou", C_NEU), ("GSE283928", C_283), ("GSE207937", C_207)]:
    m = coh == c
    ax[0, 0].scatter(PC1[m], PC2[m], s=18, c=col, alpha=0.6,
                     marker="^" if c != "GSE283928" else "o",
                     label=c, edgecolors="none")
ax[0, 0].set_xlabel("PC1 (29.2%)"); ax[0, 0].set_ylabel("PC2 (9.7%)")
ax[0, 0].set_title("(A) Joint 3-cohort PCA (687,463 common probes)")
ax[0, 0].legend(loc="upper right", fontsize=9)
ax[0, 0].grid(alpha=0.2)

# (B) projection heatmap (-distance to centroids)
def build_sim(ass):
    sim = []
    for r in ass:
        sim.append([-float(r[f"dist_c{c}"]) for c in clusters])
    sim = np.array(sim)
    order = np.argsort([int(r["assigned_cluster_eucl"]) for r in ass])
    return sim[order], order
neu_sim, neu_o = build_sim(neu_ass)
g207_sim, g207_o = build_sim(g207_ass)
sim = np.vstack([neu_sim, g207_sim])
# row color bar
rowcol = np.array([C_NEU] * len(neu_ass) + [C_207] * len(g207_ass))
im = ax[0, 1].imshow(sim, aspect="auto", cmap="viridis")
ax[0, 1].set_xticks(range(len(clusters)))
ax[0, 1].set_xticklabels([f"c{c}" for c in clusters])
ax[0, 1].set_yticks([])
ax[0, 1].set_title("(B) Similarity to GSE283928 tumor centroids\n(Neou top, GSE207937 bottom; red=NFPT-rich c0)")
# annotate cohort split
ax[0, 1].axhline(len(neu_ass) - 0.5, color="white", lw=1.5)
fig.colorbar(im, ax=ax[0, 1], fraction=0.046, pad=0.04, label="-Euclidean dist")

# (C) Moller-key methylation bar
means = [g[1].mean() for g in groups]
sds = [g[1].std(ddof=1) for g in groups]
xs = np.arange(len(groups))
ax[1, 0].bar(xs, means, yerr=sds, capsize=5, color=[C_NEU, "#ff9999", C_207, "#99cc99"],
             edgecolor="black", linewidth=0.6)
for i, (m, s) in enumerate(zip(means, sds)):
    ax[1, 0].text(i, m + s + 0.005, f"{m:.3f}", ha="center", fontsize=9)
ax[1, 0].set_xticks(xs); ax[1, 0].set_xticklabels([g[0] for g in groups], fontsize=8)
ax[1, 0].set_ylabel("Mean beta @ Moller-2025 key genes")
ax[1, 0].set_title("(C) NFPT hypermethylation replicates\n(Neou NFPT 0.687 vs Func 0.563; p<0.001)")
ax[1, 0].set_ylim(0.45, 0.75)
ax[1, 0].grid(axis="y", alpha=0.2)

# (D) lineage score boxplot
bp = ax[1, 1].boxplot([g[2] for g in groups], patch_artist=True, showmeans=True)
ax[1, 1].set_xticks(xs); ax[1, 1].set_xticklabels([g[0] for g in groups], fontsize=8)
for i, b in enumerate(bp["boxes"]):
    b.set_facecolor([C_NEU, "#ff9999", C_207, "#99cc99"][i])
    b.set_alpha(0.6)
for med in bp["medians"]:
    med.set_color("black")
ax[1, 1].set_ylabel("Lineage score  SF1 - (TPIT+PIT1)/2")
ax[1, 1].set_title("(D) Crude lineage proxy (replicates GSE207937:\nNFPT slightly less SF1-dominant)")
ax[1, 1].axhline(0, color="gray", lw=0.8, ls="--")
ax[1, 1].tick_params(axis="x", labelsize=8)
ax[1, 1].grid(axis="y", alpha=0.2)

plt.tight_layout(rect=[0, 0, 1, 0.96])
png = os.path.join(OUT, "neou_crosscohort_figure.png")
pdf = os.path.join(OUT, "neou_crosscohort_figure.pdf")
plt.savefig(png, dpi=150)
plt.savefig(pdf)
print("saved", png, pdf)
print("Moller means:", {g[0].replace("\n", " "): round(float(g[1].mean()), 3) for g in groups})
print("Lineage medians:", {g[0].replace("\n", " "): round(float(np.median(g[2])), 3) for g in groups})
