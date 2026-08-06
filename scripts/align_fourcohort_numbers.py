#!/usr/bin/env python3
# Align the four-cohort external-validation percentages in the manuscript
# to the authoritative consolidated-validation JSON
# (results/Moller2025/external_validation_consolidated/external_validation_consolidated.json):
#   Neou 2020 NFPT      26/32  = 81.2%
#   GSE207937 non-sec.  34/45  = 75.6%   <-- was wrongly 78% (35/45)
#   Moller 2025 raw     101/117 = 86.3%
#   GSE115783 NFPA       32/34  = 94.1%
#   mean(all)  = 84.3%   ; mean(EPIC-only) = 81.0%
import io, sys

path = r"C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/docs/manuscript-draft.md"
with io.open(path, "r", encoding="utf-8") as f:
    txt = f.read()

repl = {
    "mapped to GSE283928 tumor cluster 0 (Neou 81%, GSE207937 78%), whereas functioning tumors mapped to clusters 3/4 (Neou 96%).":
        "mapped to GSE283928 tumor cluster 0 (Neou 81.2%, GSE207937 75.6%), whereas functioning tumors mapped to clusters 3/4 (Neou 96%).",
    "assigned 86% to the SF1-lineage c0 cluster.":
        "assigned 86.3% to the SF1-lineage c0 cluster.",
    "94% of GSE115783 NFPA tumors mapped to c0":
        "94.1% of GSE115783 NFPA tumors mapped to c0",
    "a mean of 84.8% of NF-PitNETs/NPFAs mapped to the c0 (SF1-lineage) centroid.":
        "a mean of 84.3% of NF-PitNETs/NPFAs mapped to the c0 (SF1-lineage) centroid.",
    "independently within* Neou (81%) and GSE207937 (78%) without reference":
        "independently within* Neou (81.2%) and GSE207937 (75.6%) without reference",
    "In Neou, 26 of 32 NF-PitNETs (81%) were assigned to cluster 0":
        "In Neou, 26 of 32 NF-PitNETs (81.2%) were assigned to cluster 0",
    "In GSE207937, 35 of 45 NF-PitNETs (78%) mapped to cluster 0":
        "In GSE207937, 34 of 45 NF-PitNETs (75.6%) mapped to cluster 0",
    "the 81% assignment rate of Neou NF-PitNETs to cluster 0 was highly significant":
        "the 81.2% assignment rate of Neou NF-PitNETs to cluster 0 was highly significant",
    "the 81\u201386% c0 mapping seen across the three EPIC-850K cohorts":
        "the 75.6\u201386.3% c0 mapping seen across the three EPIC-850K cohorts",
    "101/117 (86%) of M\u00f8ller's NF-PitNETs mapped to the GSE283928 c0 (SF1-lineage) centroid":
        "101/117 (86.3%) of M\u00f8ller's NF-PitNETs mapped to the GSE283928 c0 (SF1-lineage) centroid",
    "Neou 2020 NFPT 81%, GSE207937 non-secreting 78%, M\u00f8ller raw-idat 86%) and one cross-platform HM450 cohort (GSE115783 NFPA 94%)\u2014a mean of 84.8% of NF-PitNETs/NPFAs map":
        "Neou 2020 NFPT 81.2%, GSE207937 non-secreting 75.6%, M\u00f8ller raw-idat 86.3%) and one cross-platform HM450 cohort (GSE115783 NFPA 94.1%)\u2014a mean of 84.3% of NF-PitNETs/NPFAs map",
    "three EPIC-850K: Neou 2020 NFPT 81%, GSE207937 non-secreting 78%, M\u00f8ller raw-idat 86%; one cross-platform HM450: GSE115783 NFPA 94%; overall mean 84%, no k-labels)":
        "three EPIC-850K: Neou 2020 NFPT 81.2%, GSE207937 non-secreting 75.6%, M\u00f8ller raw-idat 86.3%; one cross-platform HM450: GSE115783 NFPA 94.1%; overall mean 84.3%, no k-labels)",
    "94% of NFPA tumors (32/34) mapped to the c0 (SF1-lineage) centroid":
        "94.1% of NFPA tumors (32/34) mapped to the c0 (SF1-lineage) centroid",
    "All four converge on 78\u201394% c0 mapping (overall mean 84.8%).":
        "All four converge on 75.6\u201394.1% c0 mapping (overall mean 84.3%; 81.0% across the three EPIC-850K cohorts alone).",
}

bad = []
for old, new in repl.items():
    n = txt.count(old)
    if n != 1:
        bad.append((old, n))
        continue
    txt = txt.replace(old, new, 1)

if bad:
    sys.stderr.write("REPLACEMENTS FAILED (expected count 1 each):\n")
    for old, n in bad:
        sys.stderr.write("  count=%d  %r\n" % (n, old[:80]))
    sys.exit(1)

with io.open(path, "w", encoding="utf-8") as f:
    f.write(txt)
print("OK: %d replacements applied" % len(repl))
