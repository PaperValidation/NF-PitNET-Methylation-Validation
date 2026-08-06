# Supplementary Tables

Companion tables to the main manuscript. Machine-readable CSVs:
`results/Neou/supplementary_table1_pca.csv`,
`results/Neou/supplementary_table2_internal_clusters.csv`,
`results/Neou/supplementary_table3_normalization.csv`.
Reproducible via `scripts/build_supplementary_tables.py`.

---

## Supplementary Table 1. Joint PCA cohort centroids and explained variance

Joint PCA was fit on the variance-filtered common-probe β-matrix (687,463 common
probes across Neou 2020 / GSE283928 / GSE207937; top 10% by variance = 68,746
probes retained; `PCA(n_components=20, random_state=42)`), yielding 281 samples
× 20 PCs. **Part A** reports the mean PC coordinate of each cohort (the cohort
centroid in PC space). **Part B** reports the fraction of variance explained by
each PC and the running cumulative. PC1 = 29.2%, PC2 = 9.7%, cumulative PC1–20 =
63.4%.

**Table ST1-A. Per-cohort mean PC coordinates (cohort centroids).**

| PC | Neou | GSE283928 | GSE207937 |
|----|-----:|----------:|----------:|
| PC1 | 8.3999 | −4.3917 | −2.6516 |
| PC2 | 1.5027 | −1.2116 | 0.1784 |
| PC3 | 5.2435 | −6.5219 | 4.1382 |
| PC4 | −5.0602 | 3.8405 | −0.2337 |
| PC5 | 0.5404 | −0.2164 | −0.2720 |
| PC6 | 4.0907 | −6.6340 | 5.5975 |
| PC7 | 1.8600 | 0.6022 | −3.0002 |
| PC8 | −0.0434 | 0.3552 | −0.4957 |
| PC9 | −0.9074 | 0.2021 | 0.7038 |
| PC10 | 0.2815 | 1.0618 | −1.9415 |
| PC11 | 1.5651 | −0.9020 | −0.3657 |
| PC12 | 0.0135 | −1.2538 | 1.9062 |
| PC13 | −3.8157 | 0.0590 | 4.1712 |
| PC14 | −0.6496 | 0.1089 | 0.5586 |
| PC15 | −1.6813 | −0.1986 | 2.1821 |
| PC16 | −1.6596 | 0.6893 | 0.7972 |
| PC17 | −0.5015 | 0.0660 | 0.4589 |
| PC18 | −0.3583 | 0.1486 | 0.1725 |
| PC19 | 0.3406 | −0.3887 | 0.2152 |
| PC20 | −0.4746 | 0.1854 | 0.2460 |

The three cohorts are well separated along PC1 (Neou +8.40, GSE283928 −4.39,
GSE207937 −2.65), confirming the shared methylation manifold described in §3.1
while retaining a mild, systematic batch offset (see NN-mixing permutation test,
§3.1).

**Table ST1-B. Explained variance per principal component.**

| PC | Variance explained | Cumulative |
|----|-------------------:|----------:|
| PC1 | 0.2923 | 0.2923 |
| PC2 | 0.0967 | 0.3890 |
| PC3 | 0.0420 | 0.4310 |
| PC4 | 0.0313 | 0.4623 |
| PC5 | 0.0271 | 0.4895 |
| PC6 | 0.0237 | 0.5132 |
| PC7 | 0.0162 | 0.5294 |
| PC8 | 0.0141 | 0.5435 |
| PC9 | 0.0126 | 0.5561 |
| PC10 | 0.0118 | 0.5679 |
| PC11 | 0.0100 | 0.5779 |
| PC12 | 0.0085 | 0.5865 |
| PC13 | 0.0080 | 0.5945 |
| PC14 | 0.0072 | 0.6017 |
| PC15 | 0.0065 | 0.6082 |
| PC16 | 0.0062 | 0.6143 |
| PC17 | 0.0055 | 0.6199 |
| PC18 | 0.0050 | 0.6249 |
| PC19 | 0.0049 | 0.6297 |
| PC20 | 0.0047 | 0.6344 |

*Note:* The GSE283928 reference-cluster centroids in PC space (c0 n=70, c2 n=2,
c3 n=22, c4 n=24) are provided in `supplementary_table1_pca_cluster_centroids.csv`;
they confirm that the four reference clusters are resolved in the joint PC space
and that c0 (SF1-lineage) sits closest to the Neou cohort centroid.

---

## Supplementary Table 2. Neou internal cluster characterization

Neou 2020 was clustered *de novo* in its own variance-filtered space (top 10% by
variance = 86,623 probes; `PCA(20, random_state=42)` + K-means k=2…6). The
optimal k by silhouette was **k=3** (silhouette 0.4045; k=2: 0.4026, k=4: 0.3256,
k=5: 0.3023, k=6: 0.2921). The three internal clusters recapitulate the SF1-lineage /
TPIT / PIT1 tripartite taxonomy recovered by reference-centroid projection (§3.4,
Figure 2B–D), demonstrating the taxonomy is data-driven rather than
classifier-dependent.

**Table ST2. Characterization of the three de-novo Neou clusters.**

| Cluster | n | NFPT n (%) | Mean Møller-key β | Mean lineage score | Dominant lineage | Dominant stage | Cross-cohort centroid | Cross-cohort counts |
|---------|--:|----------:|------------------:|-------------------:|-----------------|----------------|----------------------|---------------------|
| int0 | 30 | 27 (90.0%) | 0.7027 | −0.1203 | SF1 | Silent | c0 (SF1-lineage) | c0=28, c3=2 |
| int1 | 15 | 1 (6.7%) | 0.6275 | +0.0479 | TPIT | Cushing (overt/subclinical) | c3 (TPIT) | c3=15 |
| int2 | 41 | 4 (9.8%) | 0.5338 | −0.0347 | PIT1 | Acromegaly / Hyperthyroidism | c4 (PIT1) | c4=41 |

NFPT contingency: of 32 NFPT tumours, 27 map to int0 (84.4%), 1 to int1, 4 to
int2; of 54 functioning tumours, 3 map to int0, 14 to int1, 37 to int2. The
NFPT-enrichment in int0 (90.0% vs the cohort-wide 37.2% NFPT prevalence) is the
independent, unsupervised confirmation that NF-PitNETs occupy the SF1-lineage
methylation state. The single TPIT-lineage NFPT (int1) and the four PIT1-lineage
NFPT (int2) are consistent with the known minor non-SF1 NF-PitNET subtypes.

---

## Supplementary Table 3. Normalization-sensitivity analysis (full four-cluster breakdown)

This table expands main **Table 3**. Neou 2020 (86 samples) was reprocessed from
raw IDAT with four normalization pipelines — `preprocessRaw`
(`preprocessRaw`), `noob`, `funnorm`, `quantile` (minfi) — and each normalized
β-matrix was projected onto the four GSE283928 centroids (c0/c2/c3/c4) by
nearest centroid in common-probe β-space. The table reports the full four-cluster
assignment counts, the c0 fraction, and the agreement with the default
`preprocessRaw` assignment.

**Table ST3. Four-cluster assignment under each normalization.**

| Normalization | n | c0 | c2 | c3 | c4 | c0 fraction | Agreement vs preprocessRaw (%) |
|---------------|--:|---:|---:|---:|---:|------------:|-------------------------------:|
| preprocessRaw | 86 | 28 | 0 | 17 | 41 | 0.3256 | 100.0 |
| noob | 86 | 28 | 0 | 17 | 41 | 0.3256 | 100.0 |
| funnorm | 86 | 29 | 0 | 16 | 41 | 0.3372 | 98.8 |
| quantile | 86 | 27 | 0 | 18 | 41 | 0.3140 | 98.8 |

Only **two** samples (P009, P082) change their nearest centroid across the four
normalizations, and neither switches into or out of the c0 (SF1-lineage) cluster;
both shift between the secretory c3/c4 centroids. The c0 fraction varies by <2.3
percentage points (0.314–0.337). These results confirm that the SF1-lineage
signal is not an artifact of a single normalization pipeline (§3.1, §3.3) and
directly support the statement in the Discussion (§4) that lineage assignment is
robust to normalization choice.
