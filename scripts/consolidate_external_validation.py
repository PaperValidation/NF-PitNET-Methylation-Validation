#!/usr/bin/env python3
"""Consolidate the cross-cohort external lineage (c0 / SF1-lineage) validation
that was already executed on INDEPENDENT cohorts into one publication-ready
deliverable: a unified table, a bar figure, and a JSON.

All cohorts are projected (NO k-labels) onto the SAME GSE283928 c0/c2/c3/c4
centroids. Numbers are pulled from each cohort's primary analysis_summary JSON
(authoritative), not from prose.

- Neou 2020 (E-MTAB-7762, EPIC):     NFPT subgroup -> c0
- GSE207937 (Silva-Junior 2023, EPIC): NS (=NFPT candidate) subgroup -> c0
- Møller 2025 raw idat (EPIC):         all NF-PitNET -> c0  (strict-test)
- GSE115783 (Bujko/Kober, HM450):      34 NFPA tumor -> c0  (CROSS-PLATFORM HM450)
  * 6 healthy-pituitary controls included as a built-in specificity control
    (only 33% map to c0; 4/6 map to c4).
"""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET"
OUTD = os.path.join(BASE, "results/Moller2025/external_validation_consolidated")
os.makedirs(OUTD, exist_ok=True)

neou = json.load(open(os.path.join(BASE, "results/Neou/neou_crosscohort_summary.json")))
gse  = json.load(open(os.path.join(BASE, "results/GSE207937/analysis_summary.json")))
mol  = json.load(open(os.path.join(BASE, "results/Moller2025/strict_projection/strict_projection_summary.json")))
gse115783 = json.load(open(os.path.join(BASE, "results/Moller2025/gse115783_validation/gse115783_summary.json")))

# --- Neou NFPT -> c0
n_nfpt = neou["neou_n_nfpt"]
n_c0_neou = neou["neou_nfpt_vs_cluster"]["NFPT"]["0"]
# --- GSE207937 NS -> c0
n_ns = gse["subtype_counts"]["NS"]
n_c0_gse = gse["crosscohort_subtype_vs_cluster"]["NS"]["0"]
# --- Møller all NF-PitNET -> c0
n_mol = mol["n_linked"]
n_c0_mol = mol["assigned_counts"]["c0"]
# --- GSE115783 NFPA tumor -> c0 (HM450, cross-platform)
g_tumor = gse115783["groups"]["tumor"]
n_g = g_tumor["n"]
n_c0_g = g_tumor["assignment_counts"]["c0"]

rows = [
    {"cohort": "Neou 2020", "accession": "E-MTAB-7762", "platform": "EPIC 850K",
     "subgroup": "NFPT", "n": n_nfpt, "c0": n_c0_neou,
     "c0_pct": round(100*n_c0_neou/n_nfpt, 1)},
    {"cohort": "GSE207937", "accession": "GSE207937", "platform": "EPIC 850K",
     "subgroup": "NS (non-secreting = NFPT candidate)", "n": n_ns, "c0": n_c0_gse,
     "c0_pct": round(100*n_c0_gse/n_ns, 1)},
    {"cohort": "Møller 2025 (raw idat)", "accession": "Mendeley vr2wb38r6t", "platform": "EPIC 850K",
     "subgroup": "all NF-PitNET", "n": n_mol, "c0": n_c0_mol,
     "c0_pct": round(100*n_c0_mol/n_mol, 1)},
    {"cohort": "GSE115783", "accession": "GSE115783", "platform": "HM450 (cross-platform)",
     "subgroup": "NFPA tumor", "n": n_g, "c0": n_c0_g,
     "c0_pct": round(100*n_c0_g/n_g, 1)},
]
mean_pct = round(float(np.mean([r["c0_pct"] for r in rows])), 1)
epic_rows = [r for r in rows if r["platform"].startswith("EPIC")]
epic_mean = round(float(np.mean([r["c0_pct"] for r in epic_rows])), 1)

summary = {
    "description": "Independent external validation of the c0 (SF1-lineage) methylation mapping across independent cohorts (3x EPIC-850K + 1x HM450 cross-platform), projected (no k-labels) onto the same GSE283928 c0/c2/c3/c4 centroids.",
    "centroid_source": "GSE283928 (Belakhoua 2025, EPIC 850K)",
    "rows": rows,
    "mean_c0_pct_all": mean_pct,
    "mean_c0_pct_EPIC_only": epic_mean,
    "cross_platform_HM450": {
        "cohort": "GSE115783", "platform": "Illumina HumanMethylation450",
        "tumor_n": n_g, "tumor_c0": n_c0_g, "tumor_c0_pct": round(100*n_c0_g/n_g, 1),
        "normal_n": gse115783["groups"]["normal"]["n"],
        "normal_c0": gse115783["groups"]["normal"]["assignment_counts"]["c0"],
        "normal_c0_pct": gse115783["groups"]["normal"]["c0_fraction"]*100,
        "note": "Healthy-pituitary controls map predominantly to c4 (4/6), only 33% to c0 -> confirms c0 is a TUMOR SF1-lineage identity, not a generic pituitary signature (specificity control)."
    },
    "note": "The cross-PLATFORM (HM450) confirmation via GSE115783 (34 NFPA) is now COMPLETE: 94.1% of independent NFPA tumors map to c0, replicating the EPIC-850K finding on a different array platform bridged by 333,574 common cg probes.",
}
json.dump(summary, open(os.path.join(OUTD, "external_validation_consolidated.json"), "w"), indent=2)

# --- figure (publication-ready, no overlapping labels)
fig, ax = plt.subplots(figsize=(9.5, 5.2))
# compact 2-line labels; full subgroup information lives in the title/caption
labels = [
    f"Neou 2020\nNFPT (n={rows[0]['n']})",
    f"GSE207937\nNS (n={rows[1]['n']})",
    f"Møller 2025\nall NF-PitNET (n={rows[2]['n']})",
    f"GSE115783\nNFPA tumor (n={rows[3]['n']})",
]
pcts = [r["c0_pct"] for r in rows]
colors = ["#2c7fb8", "#41b6c4", "#fec44f", "#d95f0e"]   # HM450 = orange
bars = ax.bar(labels, pcts, color=colors, edgecolor="black", linewidth=0.4)

# mean reference lines and labels (white bbox so labels stay readable over any bar)
txt_bbox = dict(boxstyle="round,pad=0.25", facecolor="white", alpha=0.92, edgecolor="none")
ax.axhline(epic_mean, ls="--", color="#2c7fb8", lw=1)
ax.text(3.55, epic_mean + 1.0, f"EPIC mean = {epic_mean}%", color="#2c7fb8", ha="right", va="bottom", fontsize=9, bbox=txt_bbox)
ax.axhline(mean_pct, ls=":", color="#666666", lw=1)
ax.text(3.55, mean_pct - 1.5, f"all-cohort mean = {mean_pct}%", color="#666666", ha="right", va="top", fontsize=9, bbox=txt_bbox)

ax.set_ylim(0, 105)
ax.set_ylabel("% of samples mapping to c0 (SF1-lineage)")
ax.set_title(
    "Independent external validation of the c0 / SF1-lineage mapping\n"
    "(3x EPIC-850K + 1x HM450 cross-platform, projected onto GSE283928 centroids, no k-labels)"
)
for b, r in zip(bars, rows):
    ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 2.0,
            f"{r['c0']}/{r['n']}\n({r['c0_pct']}%)", ha="center", va="bottom", fontsize=9)

# tighten layout and remove extra whitespace around labels
fig.subplots_adjust(left=0.10, right=0.98, bottom=0.18, top=0.86)
plt.savefig(os.path.join(OUTD, "external_validation_bar.png"), dpi=150)
plt.close()
print(json.dumps(summary, indent=2))
print(f"[done] -> {OUTD}")
