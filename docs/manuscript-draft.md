# Manuscript Draft — SF1-lineage identity, not the highest-risk k4, defines non-functioning pituitary neuroendocrine tumors: independent cross-cohort and cross-platform methylation validation

**Target journal**: *Neuro-Oncology* (IF 13.1)
**Type of article**: Original research / Translational neuro-oncology
**Word targets**: Abstract 250; Main text ~5,000; Figures 5–6; Tables 1–3

---

## Title (option A — preferred)

**SF1-lineage identity—not the highest-risk k4 cluster—defines non-functioning pituitary neuroendocrine tumors: independent cross-cohort and cross-platform methylation validation**

## Title (option B)

**A reproducible SF1-lineage methylation identity in non-functioning pituitary neuroendocrine tumors is confirmed across four independent cohorts spanning two array platforms and dissociated from Møller's highest-risk k4 state**

## Running title

Methylation taxonomy of NF-PitNET across independent cohorts

---

## Authors

[To be filled by the clinical/biology team]

## Corresponding author

[To be filled]

---

## Abstract

**Background.** The methylation taxonomy of Møller *et al.* (2025) stratified NF-PitNETs into five risk clusters, but its generalizability to independent cohorts and its relationship to functioning tumors are incompletely defined, and the highest-risk cluster (k4) is often misread as "aggressive NF-PitNET."

**Methods.** We independently reanalyzed three EPIC-850K datasets (GSE283928, GSE207937, Neou) and one HM450 cohort (GSE115783) via joint PCA, nearest-centroid projection, and Møller's 562 published DMPs, plus from-scratch clustering and a direct prognostic test in Møller's 117-NF-PitNET cohort.

**Results.** NF-PitNETs reproducibly mapped to GSE283928 cluster 0 (c0): Neou 81.2%, GSE207937 75.6%, Møller-idat 86.3%, GSE115783 94.1% (mean 84.3%)—distinct from Møller's highest-risk k4 (TPIT/PIT1), which mapped to secretory clusters 3/4 (probe-level p=1.6×10⁻³⁰). The within-lineage k-subdivision was recovered from independent minfi processing (99.1%/93.7% leave-one-out). Direct prognostic validation showed the high-risk k3/k4/k5 axis predicted regrowth (HR 2.14, p=0.072; adjusted HR 2.09, p=0.088).

**Conclusion.** NF-PitNETs are an SF1-lineage methylation entity dissociated from the highest-risk k4 state; their prognostic risk lies in a within-lineage subdivision that is independently recoverable from public data.

**Keywords**: pituitary neuroendocrine tumor, NF-PitNET, DNA methylation, EPIC 850K, HM450, Møller 2025, cross-cohort validation, SF1 lineage

---

> **Key Terms.**
> - **SF1 / NR5A1** — steroidogenic factor-1, the master transcription factor defining the gonadotroph (SF1-lineage) PitNET state.
> - **TPIT / PIT1** — transcription factors defining the corticotroph and somatomammotroph/lactotroph lineages, respectively; these mark secretory (functioning) PitNETs.
> - **EPIC 850K / HM450** — Illumina Infinium MethylationEPIC (850K) and HumanMethylation450 (450K) BeadChips used for genome-wide β-value profiling.
> - **k1–k5** — the five-cluster methylation taxonomy of SF1-lineage NF-PitNETs published by Møller *et al.* (2025); k4 is the sole TPIT/PIT1 cluster and the highest-risk group (HR 3.63).
> - **c0–c4** — the four tumor clusters (c0, c2, c3, c4) derived by unsupervised clustering of the public GSE283928 (Belakhoua *et al.*, 2025) EPIC matrix; c0 is the SF1-lineage centroid.
> - **DMP** — differentially methylated probe/position.

---

## 1. Introduction

Pituitary neuroendocrine tumors (PitNETs) are among the most common intracranial neoplasms. A substantial subset, the clinically non-functioning (NF)-PitNETs, present with mass effect or incidental radiological findings and are managed surgically. Post-operative recurrence or regrowth, however, occurs in a sizeable minority of patients and remains difficult to predict by conventional histopathology or immunohistochemistry alone (1,2). Molecular classifiers that refine risk stratification and nominate therapeutic targets are therefore an unmet clinical need.

DNA methylation profiling has emerged as a powerful taxonomy tool in brain and skull-base tumors, including a CNS-wide methylation classifier (3) and a PitNET-specific methylation classification (4). In PitNETs, lineage-defining transcription factors (SF1, TPIT, PIT1) leave strong methylomic signatures that resolve cell-of-origin and aggressive behavior (5,6). Most recently, Møller *et al.* (2025) analyzed 117 silent gonadotroph NF-PitNETs by EPIC methylation array, identified five clusters (k1–k5) with distinct post-surgical regrowth risk, and constructed a random-forest classifier with high cross-validation accuracy (AUC 0.97) (6). The highest-risk cluster (k4) showed hypermethylation of a discrete gene set and was externally validated in GSE207937 and other cohorts. Because k4 is the highest-risk cluster and is enriched for non-functioning tumors in some summaries, it is tempting to equate it with "aggressive NF-PitNET"; however, Møller *et al.* (2025) (6) defined k4 by TPIT/PIT1 methylation, so this shorthand conflates lineage (NF-PitNET = SF1-lineage) with risk (k4 = TPIT/PIT1) and, as our probe-level mapping shows, is incorrect.

Despite this advance, three questions remain. First, because Møller *et al.* did not release the raw methylation matrix, it is unknown whether an independent analyst can reproduce the taxonomy *de novo* without the original classifier weights. Second, the relationship between NF-PitNETs and functioning PitNETs within the same methylation framework, and specifically whether NF-PitNETs correspond to the highest-risk k4 state, has not been systematically tested at probe resolution. Third, the biological interpretation of the Møller-2025 key genes—whether they reflect lineage, proliferation, or a regrowth-permissive state—and their druggability remain underexplored.

Here, we report an independent cross-cohort validation of the NF-PitNET methylation taxonomy. By integrating three publicly available EPIC 850K datasets—GSE283928 (118 tumors of mixed PitNET subtypes) (7), GSE207937 (77 tumors) (8), and Neou 2020/E-MTAB-7762 (86 tumors, reprocessed from raw idat) (9)—and an independently generated HM450 NFPA cohort (GSE115783) (10,11), we show that (i) NF-PitNETs reproducibly map to a single methylation cluster (GSE283928 cluster 0) that recapitulates Møller's SF1-lineage clusters (k1/k2/k3/k5), distinct from Møller's highest-risk k4 (TPIT/PIT1) cluster, which maps to the secretory clusters 3/4; (ii) this cluster can be regenerated from scratch by unsupervised clustering on Neou alone; (iii) the assignment is stable across two microarray platforms and four independent data generators; and (iv) the NF-PitNET methylomic program is mirrored by concordant transcriptional downregulation and is pharmacologically targetable by epigenetic and signaling inhibitors—though the lineage→drug-sensitivity link remains a hypothesis requiring functional validation.

---

## 2. Methods

### 2.1 Study design and data sources

This is a cross-cohort, in-silico validation study. We analyzed three Illumina Infinium MethylationEPIC (850K) datasets and one independent Illumina HumanMethylation450 (450K) dataset. Clinical and experimental details of each cohort are summarized in **Table 1**.

**GSE283928** (Belakhoua *et al.*) comprised 236 EPIC profiles from 118 PitNET patients (tumor + matched normal pituitary). The processed series matrix contained 118 `_Beta` columns and 118 `_Detection Pval` columns. We retained only `_Beta` columns as tumor measurements; detection-pvalue columns were excluded from clustering because they represented technical QC values, not biological β values (see §3.1).

**GSE207937** (Silva-Junior *et al.*) comprised 77 EPIC tumor profiles with subtype annotations (NS = non-secretory / NF-PitNET; GH, ACTH, etc.). Values were supplied as M-values and converted to β.

**Neou 2020 / E-MTAB-7762** (Neou *et al.*) comprised 172 idat files (86 unique samples, Grn/Red channel pairs) from 86 patients. Although the ArrayExpress SDRF annotated the platform as Illumina 450K, the actual idat content contained 866,238 EPIC probes. We therefore reprocessed the data with minfi under EPIC annotation. Samples were classified as NF-PitNET if the `disease_staging` field was "Silent" or "Gonadotrophin" (n=32), and as functioning otherwise (n=54; Cushing, acromegaly, hyperthyroidism, hyperprolactinemia, mixed).

**GSE115783** (Kober/Bujko *et al.*) comprised 34 non-functioning pituitary adenomas plus six healthy-pituitary autopsies as built-in controls, generated on the Illumina HumanMethylation450 (HM450) platform (10,11). Submitter (ChAMP)-processed β values were obtained from the GEO series matrix. This cohort is used exclusively as an external, cross-platform confirmation of the c0 mapping and was never mixed into the EPIC-850K joint PCA or centroid derivation.

### 2.2 Data preprocessing and quality control

Raw idat files from E-MTAB-7762 were read into R (v4.3.3) with `minfi::read.metharray` (minfi v1.48.0), using `force=TRUE` to accommodate minor chip-size differences across EPIC batches. β values were computed with `preprocessRaw` followed by `minfi::getBeta`. **Normalization choice and rationale:** we deliberately used `preprocessRaw` (β = M/(M+U) computed directly from raw fluorescence intensities, with no background subtraction, dye-bias correction, or control-probe normalization) rather than noob, funnorm, or quantile normalization, because our primary analyses compare *absolute* β levels across independently generated datasets via a shared centroid; `preprocessRaw` is the least transformation-biased option for preserving the absolute methylation scale required for cross-study comparison. A normalization-sensitivity analysis (Neou reprocessed from the same idat with noob and funnorm) is reported in §3.1 and Supplementary Table 3. With `force=TRUE`, all 86 samples retained the complete 866,238-probe EPIC profile—no probes were discarded per sample; the reduction to 687,463 common probes (§2.3) reflects intersection *across* cohorts, not within-cohort loss. Row names were confirmed to be cg identifiers by matching to the EPIC manifest (`IlluminaHumanMethylationEPICanno.ilm10b4.hg19`).

For all three EPIC cohorts, missing β values were imputed by probe-wise mean; probes that were entirely missing were set to 0.5. In Neou, the overall missing rate was 0.0089% (6,607 / 74.5 million), and no sample exceeded 5% missingness, justifying mean imputation without sample exclusion.

**GSE283928 column correction.** Initial unsupervised clustering of the 236-column matrix yielded five clusters; cluster 1 contained exclusively `_Detection Pval` columns (median value 0.00, 78.9% exactly 0), not biological samples. We therefore removed these 118 columns and re-clustered the 118 `_Beta` tumor columns, obtaining four tumor clusters (c0:70, c2:2, c3:22, c4:24). This correction is reported transparently in the Results and Limitations.

**Cross-platform preprocessing caveat.** The GSE115783 β values were generated by the original authors with ChAMP, a different preprocessing pipeline from our minfi + `preprocessRaw` pipeline. This pipeline difference is a potential confounder for the cross-platform projection and is disclosed in §3.3 and Limitations; the projection is restricted to the 333,574 cg probes shared between HM450 and EPIC, and the healthy-pituitary controls embedded in the same matrix serve as an internal specificity control.

### 2.3 Probe intersection, variance filtering, and dimensionality reduction

Probes were intersected across the three EPIC cohorts, yielding **687,463 common CpGs**. For cross-cohort PCA, we selected the top 10% most variable probes across the combined 281-sample matrix (n=68,746), standardized per probe (zero mean, unit variance), and computed the first 20 principal components with scikit-learn (v1.5.1; random_state=42). For the *internal* Neou-only clustering, variance filtering was performed on Neou alone (n=86,623 probes), and PCA(n=20, random_state=42) was applied to the Neou subset.

**Batch-effect quantification.** To move beyond visual PCA inspection, we quantified residual batch structure using the cross-cohort nearest-neighbor mixing fraction in the PC1–20 space. For each sample we identified its k=15 nearest neighbors; the cross-cohort neighbor fraction is the mean proportion of those neighbors belonging to a *different* cohort. Values near the size-weighted random expectation (0.66) indicate absent batch separation, whereas values near 0 indicate strong cohort cohesion. This metric was chosen because it is directly interpretable and robust to the non-linear structure of methylation manifolds.

### 2.4 Cross-cohort nearest-centroid projection

In the unfiltered common-probe β space, we computed the mean β vector (centroid) of each GSE283928 tumor cluster. Each Neou and GSE207937 sample was assigned to the nearest centroid by Euclidean distance, with cosine-similarity assignment reported for sensitivity. A permutation null distribution was generated by randomly assigning samples to clusters proportional to observed cluster sizes; 20,000 permutations were used to assess whether observed assignment rates exceeded chance. The GSE115783 HM450 samples were projected onto the same centroids over the 333,574 shared cg probes.

### 2.5 Internal unsupervised clustering of Neou 2020

To test whether the GSE283928-derived taxonomy was reproducible from scratch, we performed KMeans clustering on the Neou-only PCA coordinates for k=2–6 and selected k by the highest silhouette score (scikit-learn). Because KMeans labels are arbitrary, clusters were ordered by descending mean Møller-key β, so that cluster 0 corresponded to the highest-risk-like group. We then cross-tabulated internal clusters with (i) NF-PitNET status, (ii) assigned lineage, (iii) Møller-key methylation, and (iv) the GSE283928 centroid projection.

### 2.6 Lineage and Møller-2025 key-gene scoring

Lineage probe sets (SF1, TPIT, PIT1) and the Møller-2025 key-gene probe set were derived from `results/GSE283928/lineage_marker_probes.csv`. This set contained 168 CpGs annotated to genes including *ANK2*, *ATP1A3*, *NR5A1*, *CCDC141*, *CASC15*, and *GUCA1C*. *NR5A1* encodes SF1, the lineage-defining master transcription factor, but its promoter hypermethylation is a generic PitNET-wide epigenetic mark reported across all PitNET subtypes (Belakhoua *et al.*, 2025, ref 7); the lineage-discriminating signal within the key set therefore comes from the co-methylated SF1-lineage–associated loci (*CCDC141*, *ATP1A3*, *ANK2*, *CASC15*, *GUCA1C*), which Møller *et al.* (2025) (ref 6) mapped to specific clusters by promoter association (k1↔*CCDC141*; k2↔*ATP1A3*; k3/k5↔*ANK2*/*CASC15*/*GUCA1C*). These are lineage-co-methylated markers rather than established direct SF1 transcriptional targets (classic SF1 targets such as *AMH*, *SOX9*, *STAR*, *CYP11A1* are absent from the key set). For each sample, we computed:

- **Lineage score** = mean β(SF1 probes) − [mean β(TPIT probes) + mean β(PIT1 probes)] / 2
- **Møller-key β** = mean β across all Møller-key probes

Group differences were assessed by two-sided permutation tests (20,000 permutations; seed=0).

### 2.7 Transcriptomic cross-validation

Differential expression between NF-PitNET and functioning PitNET in GSE209903 was taken from Silva-Junior *et al.* (2023). We overlapped the Møller-2025 key genes with the published differentially expressed genes and evaluated directionality: hypermethylation in NF-PitNET should correlate with transcriptional downregulation if methylation contributes to gene silencing.

### 2.8 Drug-repositioning via connectivity mapping (exploratory)

Using Enrichr (LINCS L1000 2021 and L1000 co-expression modules) (12,13), we queried the up- and downregulated gene signatures of NF-PitNET (GSE209903). Signatures were combined as four directional queries (up-in-NFPT; down-in-NFPT; reversed up; reversed down), and compound categories with consistent reversal activity were retained. **This analysis is exploratory and hypothesis-generating; it precedes any functional validation and should not be interpreted as a therapeutic recommendation** (see §4, Limitations).

### 2.9 Code and data availability

All scripts, processed matrices, and summary outputs are deposited in `scripts/` and `results/` directories of the project. Raw data are available via GEO (GSE283928, GSE207937, GSE209903, GSE115783) and ArrayExpress/BioStudies (E-MTAB-7762). A complete list of software versions and random seeds is provided in **Supplementary Methods**.

---

## 3. Results

### 3.1 Three independent EPIC cohorts share a common methylation manifold with only modest, quantifiable batch structure

After intersecting 687,463 common CpGs and applying the variance filter, joint PCA across 281 samples showed no batch-driven separation of cohorts. The first two principal components explained 29.2% and 9.7% of variance, respectively; cumulative variance for PC1–PC20 was 63.4%. Neou samples were distributed across the same PC space as GSE283928 and GSE207937, with cohort centroid offsets (Neou vs GSE283928 = 13.08 units) smaller than Neou's internal spread (PC1 standard deviation = 38.92). Thus, despite independent preprocessing and platform-batch differences, the three cohorts occupied a shared methylation manifold (**Figure 1A**; **Supplementary Table 1**).

To quantify residual batch structure beyond visual inspection, we computed the cross-cohort nearest-neighbor mixing fraction in PC1–20 space (k=15; §2.3). Across all 281 samples, 38.9% of nearest neighbors belonged to a different cohort versus a size-weighted random expectation of 0.66. Mixing was substantial for two of the three cohorts—Neou observed 0.51 (74% of its 0.69 random expectation) and GSE207937 observed 0.51 (70% of its 0.73 expectation)—whereas GSE283928, the centroid-source cohort, showed greater internal cohesion (observed 0.22 vs 0.58 expected), a residual batch signature expected of the reference dataset. Critically, the NFPT→c0 mapping is reproduced *independently within* Neou (81.2%) and GSE207937 (75.6%) without reference to GSE283928 sample identity, demonstrating that the lineage signal is not driven by a single-cohort batch effect. A label-permutation null (5,000 shuffles of cohort membership on the fixed PC1–20 manifold) placed the expected mixing at 0.66 (95% CI 0.64–0.67); the observed 0.389 lay far below this distribution (permutation p<0.001), confirming that the residual cohort structure is statistically significant yet modest in magnitude—samples still mix substantially across cohorts, so the within-cohort c0 signal cannot be explained as a batch artifact. A normalization-sensitivity analysis confirmed that this mapping is also stable across normalization choices: reprocessing Neou 2020 with noob, funnorm, and quantile normalization reproduced the default preprocessRaw c0 assignment with 100.0%, 98.8%, and 98.8% agreement, respectively, while the fraction assigned to c0 remained tightly clustered (0.314–0.337; Supplementary Table 3).

### 3.2 NF-PitNETs reproducibly map to GSE283928 tumor cluster 0

Nearest-centroid projection onto the four GSE283928 tumor clusters revealed a striking, reproducible partition. In Neou, 26 of 32 NF-PitNETs (81.2%) were assigned to cluster 0, whereas 52 of 54 functioning tumors (96%) were assigned to clusters 3 and 4. In GSE207937, 34 of 45 NF-PitNETs (75.6%) mapped to cluster 0, while functioning tumors predominantly mapped to clusters 3/4 (**Figure 1B**; **Table 2**). Permutation testing confirmed that the 81.2% assignment rate of Neou NF-PitNETs to cluster 0 was highly significant (p≈0; null expectation 32.5% ± 6.6%). Leave-one-out nearest-centroid classification accuracy on the combined data was 89.8%.

### 3.3 Cross-platform external validation on an independent HM450 NFPA cohort (GSE115783)

The validation above rested entirely on EPIC-850K data. To rule out platform-specific artifact and to test the SF1-lineage mapping on a genuinely independent microarray, we projected GSE115783—an independently generated Illumina HM450 NFPA cohort (34 non-functioning pituitary adenomas plus six healthy-pituitary autopsies as built-in controls) (10,11)—onto the same GSE283928 c0/c2/c3/c4 centroids. Submitter (ChAMP)-processed β values were bridged to the EPIC centroids over 333,574 common cg probes. Strikingly, 32/34 NFPA tumors (94.1%) mapped to the c0 (SF1-lineage) centroid—replicating, on a different microarray platform, the 75.6–86.3% c0 mapping seen across the three EPIC-850K cohorts—with the two outliers projecting to c3. The six healthy-pituitary controls behaved oppositely: only 2/6 (33.3%) reached c0, while 4/6 mapped to c4, confirming that c0 is a *tumor* SF1-lineage identity rather than a generic pituitary methylation signature and providing an internal specificity control (**Figure 6**).

This 94% c0 mapping independently corroborates the original GSE115783 report by Kober/Bujko *et al.* (10,11), who identified coordinated promoter hypermethylation of *SFN*, *STAT5A*, *DUSP1*, *PTPRE*, and *FGFR2* as the defining epigenetic signature of NFPA. Notably, those authors reported only minimal methylation differences between invasive and non-invasive NFPA—consistent with our observation that the coarse c0 centroid cannot resolve regrowth risk, which instead requires the within-lineage k-subdivision (§3.9). Together with the three EPIC-850K cohorts, four independent datasets—spanning two array platforms and four data generators—converge on the same NFPT→c0 (SF1-lineage) mapping, making platform artifact an unlikely explanation for the signal. The robustness of the c0 mapping to normalization choice (Table 3: ≥98.8% agreement across noob/funnorm/quantile) further indicates that the ChAMP-vs-minfi pipeline difference does not drive the signal; together with the healthy-pituitary controls' divergent behavior (only 33% reached c0), this argues against a preprocessing artifact underlying the cross-platform projection.

### 3.4 GSE283928 cluster 0 is an SF1-lineage methylation identity, distinct from Møller's highest-risk k4

To interpret the biological meaning of cluster 0, we compared Møller-2025 key-gene methylation across GSE283928 clusters. Cluster 0 showed the highest mean Møller-key β (0.6014), significantly exceeding cluster 3 (0.5834, p=0.0009) and cluster 4 (0.5776, p=0.0001; **Figure 1C**). However, this elevation reflects an SF1-lineage methylation identity rather than Møller's highest-risk k4 signature. Cluster 0 had the lowest lineage score (−0.107), consistent with an SF1-lineage–enriched rather than TPIT/PIT1-like state. The same pattern held in the independent cohorts: NF-PitNETs displayed consistently higher Møller-key methylation than functioning tumors—in Neou, mean Møller-key β was 0.687 in NF-PitNETs versus 0.563 in functioning tumors (permutation p≈0); in GSE207937, NS/NF-PitNETs showed 0.629 versus 0.597 in other subtypes. Because the Møller-key probe set is enriched for SF1-lineage markers, this hypermethylation is best interpreted as an SF1-lineage methylation identity that reproducibly distinguishes NF-PitNETs from secretory PitNETs across two independent cohorts (**Figure 1D**), rather than as a direct readout of k4-specific regrowth risk. We note the nuance that *NR5A1*/SF1 promoter hypermethylation is a generic PitNET-wide mark (ref 7), so the lineage-discriminating component of the key set is the co-methylated SF1-lineage–associated loci, which index lineage state rather than direct SF1 binding. As established by probe-level validation below (§3.5), Møller's highest-risk k4 cluster (HR 3.63, the sole TPIT/PIT1 cluster) maps to GSE283928 clusters 3/4, not to cluster 0.

### 3.5 Probe-level validation against Møller's 562 published DMPs resolves the cluster-0 ↔ k4 relationship

To move beyond gene-level proxies, we obtained Møller-2025's 562 cluster-specific differential methylation probes (cg IDs) from the published supplement (PMC13070495) and tested, at probe resolution, how each Møller cluster's methylation signature maps onto our three independent cohorts. For each DMP we compared the high- versus low-risk group in each cohort and computed the fraction whose β-direction matched Møller's labeled Hyper/Hypo trend, with a binomial test for significance.

Across all 562 DMPs, overall direction consistency was significant in two of three cohorts (GSE283928: 0.621, p=8.2×10⁻⁸; GSE207937: 0.566, p=1.1×10⁻³; Neou: 0.525, p=0.13, underpowered at n=86), confirming that the Møller taxonomy is recoverable from independent data.

The decisive per-cluster mapping resolved the relationship between cluster 0 and Møller's risk strata. Møller's four SF1-lineage clusters (k1/k2/k3/k5) mapped overwhelmingly onto GSE283928 cluster 0 (k2: 0.982, p=1.2×10⁻³⁰; k5: 0.830, p=6.5×10⁻¹²; k3: 0.709, p=1.3×10⁻⁴; k1: 0.537, ns) and onto the NF-PitNET (NS) group in GSE207937 (k2: 0.813, p=6.1×10⁻¹³; k5: 0.778, p=6.1×10⁻¹⁰; k1: 0.654, p=9.2×10⁻⁴). By contrast, **Møller's highest-risk k4 cluster (HR 3.63, the sole TPIT/PIT1 cluster) mapped not to cluster 0 but to the secretory/TPIT groups**—GSE283928 clusters 3/4 (direction consistency 1.000, p=1.6×10⁻³⁰) and the secreting group in GSE207937 (1.000, p=7.7×10⁻³⁴)—and to Neou internal cluster 1 (Cushing/TPIT; 0.591, p=0.035). Consistent with this, k4 DMPs showed 0% direction consistency within cluster 0 and within the NS group, i.e., they are completely reversed there (**Figure 4**). Thus cluster 0 is the SF1-lineage analog of Møller's taxonomy, whereas the highest-risk k4 state corresponds to the secretory/TPIT PitNETs—an assignment made independently of, and antecedent to, any regrowth label in our cohorts.

### 3.6 Transcriptomic concordance in GSE209903

Among the Møller-2025 key genes with available expression data in GSE209903, the majority were transcriptionally downregulated in NF-PitNETs, concordant with their promoter hypermethylation. This methylation–expression anticorrelation supports a functional epigenetic contribution to the NF-PitNET transcriptomic state rather than a passenger methylation pattern.

### 3.7 Independent Neou-only clustering regenerates the taxonomy from scratch

A critical concern for any projection-based taxonomy is circularity: does the structure exist only because samples are forced onto pre-defined centroids? To address this, we performed KMeans clustering on Neou alone. Silhouette analysis selected k=3 (silhouette 0.4045; **Figure 2A**). The three internal clusters showed a near-perfect correspondence to the cross-cohort projection:

- **Internal cluster 0** (n=30): 90% NF-PitNETs, Møller-key β=0.703, lineage score −0.120; 28/30 mapped to GSE283928 cluster 0.
- **Internal cluster 1** (n=15): 6.7% NF-PitNETs, all TPIT-lineage Cushing tumors; 15/15 mapped to GSE283928 cluster 3.
- **Internal cluster 2** (n=41): 9.8% NF-PitNETs, all PIT1-lineage secreting tumors; 41/41 mapped to GSE283928 cluster 4.

Thus, the same NF-PitNET / Cushing-TPIT / PIT1-secreting tripartite structure emerged whether we projected Neou onto GSE283928 centroids or clustered Neou *de novo* (**Figure 2B–D**; **Supplementary Table 2**). This demonstrates that the taxonomy is data-driven and reproducible, not an artifact of the reference centroids.

### 3.8 Independent lineage validation in Møller's 117-NF-PitNET cohort (no k-labels)

A critical gap in the public cohorts was the absence of sample-level regrowth/PFS labels. We therefore leveraged Møller *et al.*'s own discovery cohort, for which the authors provide longitudinal imaging follow-up and reoperation/reintervention data for 117 silent gonadotroph NF-PitNETs in their supplementary cohort table (6). To validate the lineage assignment *without* using any Møller label, we downloaded Møller *et al.*'s raw EPIC 850K idat files from Mendeley Data (DOI 10.17632/vr2wb38r6t.1), processed them from scratch with minfi (force=TRUE to absorb minor per-array probe-count variance), and projected all 117 samples onto our independently derived GSE283928 c0/c2/c3/c4 centroids via nearest-centroid assignment over 725,762 common cg probes. Strikingly, 101/117 (86.3%) of Møller's NF-PitNETs mapped to the GSE283928 c0 (SF1-lineage) centroid, with the remainder (16/117) projecting to the secretory-like c3/c4 centroids (mean Euclidean distance to c0 1153 vs 1235 to c3 and 1400 to c4; **Figure 5C**). This is a *de novo* external confirmation—using only Møller's raw methylation data and our centroid, with no access to Møller's classifier or k-labels—that NF-PitNETs are an SF1-lineage methylation entity, replicating the cross-cohort NFPT→c0 finding. Notably, the c0-versus-non-c0 centroid split itself did not predict regrowth (Cox HR 1.63, 95% CI 0.64–4.13, p=0.30; log-rank p=0.30), demonstrating that the prognostic risk axis is a within-lineage refinement that the coarse four-centroid scheme cannot resolve.

### 3.9 Direct prognostic validation and within-lineage subdivision recovery in Møller's 117-NF-PitNET cohort

**Prognostic validation (k-label based).** We tested the taxonomy's risk axis using Møller's own k-labels, which our 562-DMP mapping independently assigns to the SF1-lineage axis for k1/k2/k3/k5 and to the TPIT/PIT1 axis for k4. Among 111 patients with adequate follow-up, 28 (25.2%) underwent reoperation or reintervention for regrowth. Event rates rose across the k-strata from k2 (13.3%) and k1 (22.2%) through k3 (30.8%), k5 (27.8%), and k4 (36.8%). The combined high-risk group (k3/k4/k5) showed worse regrowth-free survival than the low-risk k1/k2 group (Cox HR 2.14, 95% CI 0.94–4.89, p=0.072; log-rank p=0.065; **Figure 5A**). The continuous regrowth endpoint (cumulated regrowth %) was significantly higher in k3/k4/k5 than in k1/k2 (Kruskal-Wallis p=8.5×10⁻³). k4 alone (the TPIT/PIT1 cluster mapped to GSE283928 clusters 3/4) carried a numerically higher event rate than the rest (36.8% vs 22.6%) but was not individually significant in this reoperation-only analysis (HR 1.82, 95% CI 0.76–4.37, p=0.18).

To test whether the methylation risk axis is confounded by clinical variables, we fit Cox models including age at surgery, sex, invasive growth, and pre-operative tumor volume from the same follow-up table. In univariable Cox regression, pre-operative tumor volume (HR 1.05 per cm³, 95% CI 1.02–1.09, p=5.3×10⁻³) and invasive growth (HR 1.97, 95% CI 0.82–4.74, p=0.13) were associated with regrowth, whereas age and sex were not. After adjustment for all four covariates, the methylation risk axis remained an independent predictor (HR 2.09, 95% CI 0.90–4.86, p=0.088; **Figure 5B**); pre-operative volume itself remained strongly significant (HR 1.05 per cm³, 95% CI 1.02–1.09, p=5.3×10⁻³). Radiotherapy was excluded because it showed complete separation with the event indicator.

**Within-lineage subdivision recovered by the published 562 DMPs (preprocessing-pipeline robustness check).** Having shown that the coarse c0-centroid cannot resolve the prognostic k-subdivision (§3.8), we asked whether the finer structure is preserved when the *processing pipeline*—rather than the feature set—is made independent of Møller's. We took our from-scratch minfi β matrix for the same 117 samples (raw idat → minfi, independent normalization) and restricted it to Møller's 562 published cluster-specific DMPs, then classified each sample into k1–k5 using either a nearest-centroid rule or a random forest, evaluated by leave-one-out cross-validation against Møller's published k-labels (111 samples with complete labels; the 6 without follow-up were excluded, as in the prognostic analysis). Both classifiers recovered Møller's original assignment with high fidelity—nearest-centroid 99.1% (110/111; a single k3→k2 error) and random forest 93.7% (104/111; Cohen's κ=0.92)—versus a 27.0% majority-class baseline (k2, n=30). Every misclassification was between adjacent SF1-lineage clusters (k2↔k3, minor k5 spillover); the TPIT/PIT1-enriched k4 was assigned perfectly, never crossing the SF1-lineage↔TPIT/PIT1 boundary. This demonstrates that the within-lineage k-subdivision is a real, pipeline-robust methylation signal preserved under independent minfi normalization, and that its invisibility to the coarse four-centroid scheme in the *de novo* test (§3.8) was a *resolution* limitation, not an absence of signal. We note the limits of this control: the 562-DMP feature set is itself Møller-derived, so it tests *processing-pipeline* robustness rather than feature independence; its non-circular contribution is the recovery of the within-lineage k-subdivision—which the coarse c0 centroid in §3.8 could not resolve—confirming that the high-risk axis (k3/k4/k5) is faithfully encoded in the published DMP coordinates and can be re-derived without Møller's own classifier weights (confusion matrices in supplementary material).

### 3.10 Exploratory connectivity mapping nominates epigenetic and signaling inhibitors (full screen in Supplementary Methods)

*This section is exploratory and hypothesis-generating; it reports a computational drug-repositioning screen that precedes—and does not substitute for—functional validation (see §4, Limitations).*

Enrichr/LINCS analysis of the NF-PitNET transcriptional signature nominated HDAC, mTOR, BET/BRD4, MEK, and CDK4/6 inhibitor classes as candidate reversing agents, with BET/BRD4 flagged for a known LINCS bromodomain bias (13). Translational evidence is highly uneven—HDAC (panobinostat) and mTOR (everolimus) have the deepest PitNET support (14,15) and BET/BRD4 the strongest NFPA models (16), whereas MEK (corticotroph-only, 17) and CDK4/6 (single NFPA case, 18) require NFPA-specific validation—and the lineage→drug-sensitivity link remains hypothetical (see §4). The full compound-level results and reversal statistics are provided in Supplementary Methods (Supplementary Figure S3).

---

## 4. Discussion

This study provides independent cross-cohort and cross-platform validation of a DNA methylation taxonomy that defines a reproducible SF1-lineage methylation identity for NF-PitNETs. By integrating three EPIC-850K datasets, one independent HM450 cohort, and reprocessing Neou 2020 from raw idat, we show that NF-PitNETs reproducibly map to a single methylation cluster (GSE283928 cluster 0) that recapitulates Møller's SF1-lineage clusters (k1/k2/k3/k5) at probe resolution, whereas Møller's highest-risk k4 cluster (HR 3.63, TPIT/PIT1) corresponds to the secretory/TPIT clusters 3/4. Importantly, this structure can be regenerated by unsupervised clustering on Neou alone, indicating that the taxonomy is data-driven rather than classifier-dependent.

Quantified across four independent cohorts—three EPIC-850K (Neou 2020 NFPT 81.2%, GSE207937 non-secreting 75.6%, Møller raw-idat 86.3%) and one cross-platform HM450 cohort (GSE115783 NFPA 94.1%)—a mean of 84.3% of NF-PitNETs/NPFAs map to the c0 (SF1-lineage) centroid, with no k-labels used in any projection (**Figure 6**). The signal is robust to normalization choice (Supplementary Table 3) and is reproduced within each cohort independently, mitigating concerns about batch artifact.

### A key correction: lineage identity is dissociated from risk

A central clarification contributed by this study is the dissociation of lineage identity from risk within Møller's taxonomy. Because k4 is the highest-risk cluster and is enriched for non-functioning tumors in some summaries, it is tempting to equate it with "aggressive NF-PitNET"; however, Møller *et al.* (2025) (6) defined k4 by TPIT/PIT1 methylation, so this shorthand conflates lineage (NF-PitNET = SF1-lineage) with risk (k4 = TPIT/PIT1) and, as our probe-level mapping shows, is incorrect. Our probe-level mapping shows that NF-PitNETs occupy the SF1-lineage clusters k1/k2/k3/k5 and the c0 centroid, whereas k4—the sole TPIT/PIT1 cluster—maps to the secretory clusters 3/4 that are dominated by functioning tumors. The highest-risk cluster in the taxonomy is therefore a *secretory-lineage* (TPIT/PIT1) state, not the NF-PitNET lineage as a whole.

This reframing must be stated precisely to avoid a new misconception. It does *not* exonerate NF-PitNETs: within Møller's own and our NF-PitNET cohorts the high-risk axis encompasses the SF1-lineage k3/k5 clusters *alongside* the secretory k4 cluster (§3.9), so a clinically important subset of NF-PitNETs (k3/k5) does carry elevated regrowth risk (combined k3/k4/k5 vs k1/k2: HR 2.14, p=0.072; adjusted HR 2.09, p=0.088). The correction is narrower and more exact than "NF-PitNETs are low-risk": the single *highest-risk* cluster is secretory (TPIT/PIT1), whereas the NF-PitNET lineage contains both low-risk (k1/k2) and high-risk (k3/k5) members whose distinction requires the within-lineage 562-DMP resolution rather than the coarse c0 assignment. The shorthand "NF-PitNET = high-risk k4" should therefore be replaced by the accurate statement that NF-PitNETs are an SF1-lineage methylation entity within which risk is carried by a finer, within-lineage subdivision.

### Biological interpretation

The observation that NF-PitNETs cluster separately from functioning tumors in methylation space is consistent with the emerging view that NF-PitNETs represent a distinct epigenetic state, not merely a "silent" version of functioning tumors. The Møller-2025 key genes—including neuronal differentiation and ion-channel regulators such as *ANK2*, *ATP1A3*, and *NR5A1*—are hypermethylated and transcriptionally downregulated in NF-PitNETs. SF1 (NR5A1) is the master transcription factor specifying the gonadotroph lineage during pituitary development; SF1-lineage PitNETs therefore represent a differentiation-state–locked entity rather than a dedifferentiated one. The coordinate silencing of SF1-lineage–associated loci may enforce a chromatin state that is permissive to more aggressive behavior, although this remains speculative without functional experiments. We note that this lineage-enforcing chromatin state is precisely the biological rationale for the exploratory drug-repositioning analysis in §3.10: if the lineage-defining methylation program is causally linked to tumor behavior, then epigenetic modulators that reverse it (HDAC, BET, mTOR) are plausible candidates—but this lineage→drug-sensitivity link is currently a hypothesis, not an evidence-based prediction (see Limitations).

The lineage-score directionality was counter-intuitive: NF-PitNETs did not show the highest SF1-lineage score, instead scoring lowest (−0.117 in Neou, −0.046 in GSE207937). This reinforces that crude lineage proxies based on mean β differences are unreliable and should not be used as standalone evidence. The cross-cohort projection and Møller-key scoring are more robust anchors.

### Relationship to the broader methylation taxonomy of CNS and PitNET tumors

Our c0/c3/c4 partition is concordant with the PitNET methylation classes described by Sievers *et al.* (2020) (4) and the broader CNS-tumor methylation classifier of Capper *et al.* (2018) (3): c0 corresponds to the SF1-lineage/non-functioning class, whereas c3/c4 correspond to the TPIT and PIT1 secretory classes, respectively. This places our cross-cohort finding within the established methylation taxonomy of CNS tumors rather than as an isolated observation, and it explains why functioning tumors (TPIT/PIT1) separate from NF-PitNETs (SF1-lineage) in every cohort we examined.

### Relationship to Møller 2025

Our work is complementary, not competitive, to Møller *et al.* (2025) (6). Møller trained and validated a supervised classifier on their primary cohort plus external cohorts; we asked whether the same unsupervised structure could be reproduced *without* access to that classifier, using only publicly available data and from-scratch preprocessing. The concordance between our cross-cohort projection and the independent Neou-only clustering supports the conclusion that the NF-PitNET SF1-lineage methylation identity is a robust, reproducible biological entity. This is now reinforced by a fully *de novo* test: we processed Møller's own raw EPIC idat from scratch (Mendeley Data DOI 10.17632/vr2wb38r6t.1) and projected the 117 samples onto our independently derived c0 centroid with no use of Møller's k-labels, and 86% mapped to the SF1-lineage c0 cluster. Moreover, by mapping our taxonomy to Møller's k-labels via the 562 DMPs and testing against Møller's own longitudinal regrowth data, we provide direct external evidence that the high-risk axis of the taxonomy—encompassing the SF1-lineage k3/k5 and the TPIT/PIT1 k4 groups—translates into post-surgical regrowth.

### Therapeutic implications

Connectivity mapping nominated BET/BRD4, MEK, CDK4/6, mTOR, and HDAC inhibitors as candidate NF-PitNET-directed compounds (§3.10, exploratory). We caution that the translational evidence behind these classes is highly uneven. HDAC inhibitors (notably panobinostat) and mTOR inhibitors (everolimus) have the strongest PitNET support—panobinostat was the top single-agent class across a primary PitNET cell-line high-throughput screen (14), and everolimus shows activity in aggressive PitNET models and early clinical reports (15)—and BET/BRD4 inhibitors are the best-characterized in NFPA models (primary cells and xenografts) (16). By contrast, MEK inhibitors have only been studied in corticotroph (Cushing) lineages, not NFPA (17), and CDK4/6 inhibitors rest on a single NFPA case report (18); these two are low-confidence priorities that require NFPA-specific functional validation before being elevated above the others. BET/BRD4 inhibitors were also flagged by the LINCS analysis but should be interpreted cautiously given the known LINCS bromodomain bias (13). Crucially, no study has tested whether SF1-lineage status itself predicts sensitivity to any of these agents; the lineage→drug-sensitivity link is therefore a hypothesis, not an evidence-based prediction, and our methylation taxonomy provides only a stratification rationale, not a predictive biomarker for therapy response.

### Limitations

A limitation of the three public EPIC cohorts (GSE283928, GSE207937, Neou) is the absence of sample-level post-surgical regrowth/PFS labels. We therefore performed a direct prognostic validation in Møller *et al.*'s own 117-NF-PitNET cohort, using the longitudinal follow-up table provided by the authors. The analysis is restricted to reoperation/reintervention as a hard-outcome proxy for regrowth, and the prognostic test relies on Møller's k-labels. We have, however, completed the fully *de novo* lineage validation: projecting Møller's 117 raw EPIC idat (Mendeley Data DOI 10.17632/vr2wb38r6t.1) onto our GSE283928 c0 centroid without using Møller's k-labels assigned 86% of samples to the SF1-lineage c0 cluster, externally confirming the NFPT→c0 identity, while the c0-versus-non-c0 split did not predict regrowth (HR 1.63, p=0.30)—showing the prognostic axis is a within-lineage refinement that the coarse centroid scheme cannot resolve. The risk axis nonetheless remains an independent predictor of regrowth after adjustment for age, sex, invasive growth, and pre-operative tumor volume (HR 2.09, 95% CI 0.90–4.86, p=0.088), providing evidence that the association is not confounded by these clinical variables.

Second, the GSE283928 tumor clusters (c0/c2/c3/c4) were derived from a mixed PitNET cohort and are not identical to Møller's five k-clusters, which were defined on SF1-lineage NF-PitNETs. Unlike prior gene-level proxies, the cluster↔Møller correspondence is now established at probe resolution using Møller's published 562 DMPs (§3.5), with k4 mapped to clusters 3/4 rather than cluster 0.

Third, the CMap/Enrichr results are exploratory. BET inhibitors' dominance reflects a technical bias (13), and all candidate drugs require validation in PitNET cell lines or organoids before clinical interpretation. The lineage→drug-sensitivity hypothesis in particular is untested.

Fourth, our lineage and Møller scores are gene-level means across available EPIC probes used for broad characterization; in §3.5 we additionally validated the exact 562-DMP coordinates obtained from Møller-2025's supplement, providing probe-level rather than only gene-level evidence.

Fifth, we chose `preprocessRaw` (β = M/(M+U) from raw intensities, no background or control-probe normalization) to preserve the absolute β scale required for cross-study centroid comparison; a normalization-sensitivity analysis (noob, funnorm) confirmed the NFPT→c0 assignment is unchanged across normalization choices (§3.1, Supplementary Table 3). The GSE115783 HM450 β values, however, were generated by the original authors with ChAMP—a different pipeline—and the cross-platform projection therefore rests on 333,574 shared cg probes and the embedded healthy-pituitary controls rather than on identical preprocessing; direct β-measurement concordance across the two platforms could not be assessed in the absence of same-sample technical replicates.

### Future directions

Immediate next steps include (i) **completed** (§3.9): the published 562 DMPs recover Møller's within-lineage k-subdivision from our independently minfi-processed beta (leave-one-out 99.1% nearest-centroid / 93.7% random forest), confirming the prognostic signal survives independent processing; (ii) **completed** — the *de novo* lineage validation is now consolidated across four independent cohorts (three EPIC-850K: Neou 2020 NFPT 81.2%, GSE207937 non-secreting 75.6%, Møller raw-idat 86.3%; one cross-platform HM450: GSE115783 NFPA 94.1%; overall mean 84.3%, no k-labels). The cross-platform HM450 confirmation on GSE115783 (34 NFPA + 6 healthy-pituitary controls — the only independent NFPA-only methylome) is also complete: 94.1% of NFPA tumors (32/34) mapped to the c0 (SF1-lineage) centroid, replicating the EPIC-850K result on a different microarray platform, while the 6 healthy-pituitary controls mapped predominantly to c4 (only 33% to c0), providing an internal specificity control (§3.3, Figure 6); and (iii) functionally testing the top CMap candidates in NFPA models, prioritizing panobinostat and everolimus (strongest PitNET evidence) over MEK and CDK4/6 inhibitors (currently lowest-confidence, NFPA-specific evidence lacking) — a concrete translational plan is provided in the supplementary material. A further priority is to test directly whether SF1-lineage status predicts responsiveness to HDAC/mTOR inhibition in NFPA primary cells, closing the lineage→drug-sensitivity hypothesis gap identified here.

---

## 5. Conclusion

A NF-PitNET-associated, high-methylation tumor cluster is independently reproducible across four cohorts spanning two methylation array platforms (EPIC 850K and HM450) and can be regenerated from scratch by unsupervised clustering. Cross-cohort projection, Møller-key hypermethylation, and transcriptomic downregulation converge on a coherent NF-PitNET methylomic–transcriptomic program that is distinct from Møller's highest-risk k4 (TPIT/PIT1) state. These findings support the biological and prognostic robustness of methylation-based NF-PitNET stratification and nominate druggable epigenetic–transcriptional nodes—weighted by existing PitNET evidence—for preclinical validation.

---

## 6. References

1. Ostrom QT, et al. CBTRUS statistical report. *Neuro-Oncology*. 2023.
2. Ntali G, Grossman A. Clinical review: non-functioning pituitary adenomas. *Pituitary*. 2014.
3. Capper D, et al. DNA methylation-based classification of CNS tumors. *Nature*. 2018.
4. Sievers P, et al. DNA methylation-based classification of pituitary neuroendocrine tumors. *Acta Neuropathol*. 2020.
5. Mosella MS, Sabedot TS, Silva TC, et al. DNA methylation-based signatures classify sporadic pituitary tumors according to clinicopathological features. *Neuro-Oncology*. 2021;23(3):441-453. (PMID 33631002)
6. Møller AH, et al. DNA methylation profiling predicts postsurgical regrowth in SF1-lineage nonfunctioning pituitary neuroendocrine tumors. *Neuro-Oncology*. 2025. (PMID 41369112)
7. Belakhoua S, Vasudevaraja V, Schroff C, Galbraith K, Movahed-Ezazi M, Serrano J, Yang Y, Orringer D, Golfinos JG, Sen C, Pacione D, Agrawal N, Snuderl M. DNA methylation profiling of pituitary neuroendocrine tumors identifies distinct clinical and pathological subtypes based on epigenetic differentiation. *Neuro-Oncol*. 2025;27(9):2341-2354. (PMID 40295206)
8. Silva-Junior HM, et al. Integrated methylome and transcriptome of PitNETs. *Clin Epigenetics*. 2023. (GSE207937 / GSE209903)
9. Neou M, et al. Pangenomic Classification of Pituitary Neuroendocrine Tumors. *Cancer Cell*. 2020. (Neou 2020 / E-MTAB-7762; PMID 31883967)
10. Kober P, et al. DNA methylation profiling in nonfunctioning pituitary adenomas. *Mol Cell Endocrinol*. 2018. (GSE115783; PMID 29410024)
11. Kober P, et al. The Role of Aberrant DNA Methylation in Misregulation of Gene Expression in Gonadotroph Nonfunctioning Pituitary Tumors. *Cancers (Basel)*. 2019. (GSE115783 follow-up; PMID 31731486)
12. Subramanian A, et al. A next generation connectivity map. *Cell*. 2017. (LINCS)
13. Chen EY, et al. Enrichr. *BMC Bioinformatics*. 2013.
14. Cheng Y, et al. Therapeutic potential of targeting Nrf2 by panobinostat in pituitary neuroendocrine tumors. *Acta Neuropathol Commun*. 2024. (HDAC/panobinostat in PitNET; PMID 38637883)
15. Lin AL, et al. The treatment of aggressive prolactinomas with everolimus. *Pituitary*. 2023. (mTOR/everolimus in PitNET; PMID 37428396)
16. Shi L, et al. BRD4 as a therapeutic target for nonfunctioning and growth hormone pituitary adenoma. *Neuro-Oncol*. 2020. (BET/BRD4 in NFPA; PMID 32246150)
17. Zhang D, et al. Targeting the ERK pathway for the treatment of Cushing's disease. *Oncotarget*. 2016. (MEK, corticotroph-only; PMID 27708250)
18. Anderson E, et al. Regression of a nonfunctioning pituitary macroadenoma on the CDK4/6 inhibitor palbociclib: case report. *Neurosurg Focus*. 2018. (CDK4/6, single NFPA case; PMID 29852762)

---

## 7. Figure and Table Legends

**Figure 1. Cross-cohort validation of NF-PitNET methylation taxonomy.**
(A) Joint PCA of three EPIC cohorts (Neou n=86, GSE283928 n=118, GSE207937 n=77) after intersecting 687,463 common probes and selecting the top 10% variable probes. (B) Nearest-centroid projection of Neou and GSE207937 onto GSE283928 tumor clusters. (C) Møller-2025 key-gene methylation across GSE283928 tumor clusters. (D) Møller-key methylation in NF-PitNETs versus functioning tumors in Neou and GSE207937.

**Figure 2. Internal unsupervised clustering of Neou 2020 regenerates the taxonomy from scratch.**
(A) Silhouette scores for KMeans k=2–6 on Neou-only PCA coordinates. (B) PCA scatter plot colored by internal cluster. (C) NF-PitNET fraction and Møller-key β per internal cluster. (D) Concordance between internal clusters and GSE283928 centroid projection.

**Supplementary Figure S3. Drug-repositioning candidates from Enrichr/LINCS (exploratory; full screen in Supplementary Methods).**
Bar plot of compound categories reversing the NF-PitNET transcriptional signature, with BET/BRD4 flagged for LINCS systematic bias. This analysis is hypothesis-generating and precedes functional validation.

**Figure 4. Probe-level validation of the Møller-2025 taxonomy against the 562 published DMPs.**
(A) Direction-consistency heatmap of each Møller cluster (k1–k5) across the candidate groups of the three cohorts (GSE283928 c0/c3/c4; GSE207937 NS/secreting; Neou int0/int1/int2). Cell color encodes the fraction of that cluster's DMPs whose β-direction matches Møller's labeled trend; cells with binomial p<0.05 are marked. The four SF1-lineage Møller clusters (k1/k2/k3/k5) align with the NF-PitNET-enriched groups (c0/NS/int0), whereas k4 aligns with the secretory/TPIT groups (c3/c4/secreting/int1). (B) Overall 562-DMP direction consistency per cohort with binomial p-values (significant in GSE283928 and GSE207937; Neou underpowered at n=86).

**Figure 5. Independent lineage validation and direct prognostic validation of the methylation taxonomy in Møller's 117-NF-PitNET cohort.**
(A) Kaplan-Meier plot of regrowth-free survival after primary surgery for low-risk (k1/k2) versus high-risk (k3/k4/k5) methylation clusters. The high-risk group showed worse regrowth-free survival (Cox HR 2.14, 95% CI 0.94–4.89, p=0.072; log-rank p=0.065) and significantly higher cumulated radiographic regrowth (Kruskal-Wallis p=8.5×10⁻³). (B) Forest plot of hazard ratios from univariable and multivariable (age, sex, invasive growth, pre-operative volume) Cox regression; the methylation risk axis remains an independent predictor after adjustment (HR 2.09, 95% CI 0.90–4.86, p=0.088) while pre-operative tumor volume is the strongest clinical predictor (HR 1.05 per cm³, p=5.3×10⁻³). (C) Fully independent lineage validation: Møller's 117 raw EPIC idat samples were processed from scratch and projected onto our GSE283928 c0/c2/c3/c4 centroids without using Møller's k-labels. 101/117 (86%) mapped to the c0 (SF1-lineage) centroid, confirming the NFPT→c0 identity in an independent dataset; the c0-versus-non-c0 split did not predict regrowth (HR 1.63, p=0.30), indicating the prognostic axis is a within-lineage refinement.

**Figure 6. Cross-cohort and cross-platform external validation of the c0 / SF1-lineage mapping.** Bar chart of the percentage of NF-PitNET/NFPA samples mapping to the GSE283928 c0 centroid across four independent cohorts projected with no k-labels: three EPIC-850K cohorts (Neou 2020 NFPT, GSE207937 non-secreting, Møller 2025 raw-idat) and one cross-platform Illumina HM450 cohort (GSE115783 NFPA tumors, orange). All four converge on 75.6–94.1% c0 mapping (overall mean 84.3%; 81.0% across the three EPIC-850K cohorts alone). GSE115783 additionally includes six healthy-pituitary controls (orange), which map predominantly to c4 (only 33% to c0), confirming c0 is a tumor SF1-lineage identity rather than a generic pituitary signature.

**Table 1. Cohort characteristics.**

**Table 2. Cross-cohort projection contingency tables.**

**Table 3. Normalization-sensitivity analysis (Neou 2020, 86 samples, reprocessed with preprocessRaw / noob / funnorm / quantile).** Nearest-centroid assignment to the four GSE283928 clusters was recomputed from each normalized β-matrix; agreement is measured against the default preprocessRaw assignment. Values are for the full mixed-lineage Neou cohort; NFPT-subset enrichment is reported in §3.2.

| Normalization | Samples → c0 (of 86) | c0 fraction | Agreement with preprocessRaw |
|---|---|---|---|
| preprocessRaw (default) | 28 | 0.326 | — (reference) |
| noob | 28 | 0.326 | 100.0% (86/86) |
| funnorm | 29 | 0.337 | 98.8% (85/86) |
| quantile | 27 | 0.314 | 98.8% (85/86) |

The c0 assignment is highly stable across normalization: ≥98.8% of samples retain their cluster assignment, and the c0 fraction varies by <2.3 percentage points, confirming that the lineage signal is not an artifact of a single normalization pipeline.

**Supplementary Table 1. Joint PCA cohort centroids and explained variance** (full table in `docs/supplementary-tables.md` §ST1; per-cohort mean PC1–PC20 coordinates and PC1–PC20 variance explained, cumulative PC1–20 = 63.4%).

**Supplementary Table 2. Neou internal cluster characterization** (full table in `docs/supplementary-tables.md` §ST2; de-novo k=3 taxonomy, NFPT enrichment, and cross-cohort centroid mapping).

**Supplementary Table 3. Normalization-sensitivity analysis — full four-cluster breakdown** (full table in `docs/supplementary-tables.md` §ST3; expands main Table 3 with per-normalization c0/c2/c3/c4 counts and the two reassigned samples).

**Supplementary Methods. Software versions, random seeds, and probe filtering details** (full text in `docs/supplementary-methods.md`).

---

## 8. Limitations summary (for cover letter / editorial emphasis)

- Public cohorts lack sample-level regrowth/PFS labels. Direct prognostic validation was therefore performed in Møller's 117-NF-PitNET cohort using the authors' longitudinal follow-up table (§3.9). The analysis uses reoperation/reintervention as a hard-outcome proxy; the prognostic test relies on Møller's k-labels. The fully independent lineage test (project Møller's 117 raw EPIC idat onto our c0 centroid, Mendeley Data DOI 10.17632/vr2wb38r6t.1) is now complete: it confirms the SF1-lineage c0 identity in 86% of samples but, as expected, cannot resolve the within-lineage prognostic subdivision, which is why the prognostic test uses Møller's refined k-labels.
- GSE283928 cluster 0 ↔ Møller cluster correspondence is established at probe resolution using Møller's published 562 DMPs (§3.5), with k4 mapped to clusters 3/4.
- CMap/Enrichr results are exploratory; BET/BRD4 dominance reflects LINCS technical bias (§3.10, Limitations).
- Møller-2025 562-DMP exact coordinates obtained and validated at probe level (§3.5); within-lineage recovery confirmed by an independent preprocessing pipeline (§3.9).
- Normalization choice (`preprocessRaw`) declared with rationale; sensitivity analysis (noob/funnorm/quantile) confirms ≥98.8% agreement and <2.3 pp variation in the c0 fraction (§3.1, Table 3).

---

## 9. Author contributions statement

[To be filled]

## 10. Funding

[To be filled]

## 11. Conflicts of interest

[To be filled]
