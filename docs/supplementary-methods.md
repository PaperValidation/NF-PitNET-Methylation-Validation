# Supplementary Methods

This document provides (S1) the software versions, random seeds, and probe-filtering
details referenced in the main text, and (S2) the full-screen drug-repositioning
analysis that §3.10 summarizes. The companion figure is **Supplementary Figure S3**
(`results/CMap/supplementary_figure_s3.png`).

---

## S1. Software versions, random seeds, and probe-filtering details

### S1.1 Software stack

| Component | Version | Role |
|---|---|---|
| R | 4.3.3 | Raw `.idat` processing (Mendeley re-processing) |
| Bioconductor `minfi` | 1.48.0 | Illumina EPIC/HM450 β-value extraction (`read.metharray`, `force=TRUE`) |
| `IlluminaHumanMethylationEPICmanifest` / `EPICanno.ilm10b4.hg19` | Bioconductor release aligned to R 4.3 | EPIC-850K probe annotation |
| `IlluminaHumanMethylation450kmanifest` | Bioconductor release aligned to R 4.3 | HM450 probe annotation (GSE115783 bridging) |
| Python | 3.13 | All downstream statistics, PCA, projection, figures |
| `scikit-learn` | 1.5.1 | PCA, `NearestNeighbors` (k=15 mixing), KMeans, `RandomForestClassifier` (LOO) |
| `pandas` / `numpy` | 2.x / 1.26.x | Tabular I/O and numerics |
| `matplotlib` / `scipy` | 3.8.x / 1.11.x | Plotting and permutation statistics |
| `gseapy` | 1.3.1 | Enrichr / LINCS-L1000 API client (connectivity mapping) |

### S1.2 Random seeds and deterministic settings

- **Joint / Neou-only PCA:** `random_state = 42` (scikit-learn `PCA`, 20 components).
- **Permutation nulls (§2.4, §2.6):** 20,000 label shuffles, `seed = 0`.
- **Nearest-neighbor mixing (§2.3):** `k = 15`; deterministic under fixed PC1–20 coordinates.
- **Internal KMeans (§2.5):** `k = 2..6`, cluster selection by maximum silhouette (selected **k=3**, silhouette 0.4045); labels ordered by descending mean Møller-key β.
- **Leave-one-out 562-DMP classification (§3.9):** deterministic; 111 labelled samples, 6 unlabelled (no follow-up) excluded.

### S1.3 Probe intersection and variance filtering (§2.3)

- Probes were intersected across the three EPIC-850K cohorts (GSE283928, GSE207937, Neou), yielding **687,463 common CpGs** after manifest-matched QC (detection-p filter, removal of XY, polymorphic, and multi-mapping probes per standard `minfi`/`ChAMP` QC).
- **Cross-cohort PCA:** the **top 10% most variable probes** across the combined 281-sample matrix (n = 68,746) were retained, standardized per probe (zero mean, unit variance), and reduced to the first 20 principal components.
- **Neou-only internal clustering:** variance filtering performed on Neou alone (n = 86,623 probes); PCA(n=20) applied to the Neou subset.
- **GSE115783 (HM450) bridging (§3.3):** projected onto GSE283928 centroids over **333,574 shared cg probes** (submitter ChAMP-processed β values).

---

## S2. §3.10 full screen — Exploratory connectivity mapping (Enrichr/LINCS)

> This analysis is **exploratory and hypothesis-generating**. It reports a computational
> drug-repositioning screen that *precedes* — and does not substitute for — functional
> validation. The lineage→drug-sensitivity link remains a hypothesis (see main text §4 and Limitations).

### S2.1 Disease signature source

- Derived from `GSE209903` (FPT_vs_NFPT DESeq2 output: `log2FoldChange` / `padj`).
- Contrast is `FPT_vs_NFPT`, so **NFPT-up = log2FC < 0** (genes *higher* in NFPT).
- Lenient set (padj < 0.05): NFPT-up **4,717** / NFPT-down **5,135** Ensembl genes.
- Stringent set (padj < 0.05 and |log2FC| ≥ 1): NFPT-up **1,552** / NFPT-down **1,917** Ensembl genes (used for enrichment).
- Ensembl→symbol mapping via `mygene` (human): stringent set → **1,431 up / 1,771 down** symbols.

### S2.2 Connectivity framework

- The classic CMap `build02` rank-matrix (~6,100 instances / 1,309 compounds) was unreachable
  in the analysis environment (Google Storage 404; GenePattern / `cmapR` endpoints failing).
- We therefore adopted **Enrichr's LINCS/L1000 drug-perturbation gene sets** as an equivalent
  connectivity framework — the same "disease signature vs drug-perturbation signature" question,
  and a standard current drug-repositioning approach, with no ~1 GB download required.
- Tool: `gseapy 1.3.1` + Enrichr API (routed via SOCKS5 proxy).

### S2.3 Reversal score definition

For each LINCS directional perturbation library we computed four enrichments and combined them
with the Enrichr Combined Score:

| Query × library | Biological meaning | Sign |
|---|---|---|
| NFPT-up × `Chem_Pert_down` | drug *down-regulates* disease-up genes | reversal (+) |
| NFPT-down × `Chem_Pert_up` | drug *up-regulates* disease-down genes | reversal (+) |
| NFPT-up × `Chem_Pert_up` | drug further up-regulates disease-up genes | pro-disease (−) |
| NFPT-down × `Chem_Pert_down` | drug further down-regulates disease-down genes | pro-disease (−) |

`reversal = (up↓ + down↑) − (up↑ + down↓)`; higher = more likely to reverse the NFPT program.
A direction-agnostic overlap library (`Drug_Perturbations_from_GEO_2014`, classic CMap style) was
used as a supplementary check.

### S2.4 Compound-level results

**Drug-class frequency among all reversal-positive perturbations** (n = **15,452** LINCS instances;
"Other / non-prioritized" = 11,283 instances excluded from the ranked list below):

| Drug class | Reversal-positive instances |
|---|---:|
| BET/BRD (bromodomain) | 2,713 |
| HSP90 | 359 |
| HDAC | 288 |
| SFK/ABL | 239 |
| mTOR/PI3K | 223 |
| EGFR/HER | 118 |
| CDK | 101 |
| MEK | 61 |
| RAF/MEK | 32 |
| Glucocorticoid | 18 |
| Hormone | 9 |
| JAK/STAT | 8 |

**Top 20 individual reversal candidates** (by reversal score; compound class annotated by name matching):

| # | Compound (cleaned) | Class | Reversal score |
|---:|---|---|---:|
| 1 | BRD-K45205755 | BET/BRD | 48,972.5 |
| 2 | gr-235 | Other | 48,972.4 |
| 3 | pirenperone | Other | 254.5 |
| 4 | WZ-4-145 | Other | 170.4 |
| 5 | dasatinib | SFK/ABL | 132.0 |
| 6 | BRD-K21170640 | BET/BRD | 112.5 |
| 7 | BRD-K27170250 | BET/BRD | 110.8 |
| 8 | salmeterol | Other | 110.8 |
| 9 | s1173 | Other | 109.2 |
| 10 | PD-0325901 | MEK | 105.1 |
| 11 | selumetinib | MEK | 93.8 |
| 12 | AZD-8055 | mTOR/PI3K | 90.4 |
| 13 | AZD-7762 | Other | 88.6 |
| 14 | foretinib | EGFR/HER | 88.6 |
| 15 | chloride | Other | 86.7 |
| 16 | CGP-60474 | Other | 83.1 |
| 17 | nilotinib | SFK/ABL | 83.1 |
| 18 | JW-7-24-1 | Other | 80.6 |
| 19 | vemurafenib | RAF/MEK | 79.5 |
| 20 | ag | Other | 77.6 |

*Note:* the two dominant instances (BRD-K45205755, gr-235) carry extreme Combined Scores and
reflect the known LINCS bromodomain bias; conclusions are drawn at the **drug-class** level, not
from single instances. Compounds labelled "Other" are unmapped LINCS/ChEMBL identifiers.

### S2.5 Evidence-tier annotation (PitNET translational evidence)

The screen nominates classes; translational evidence in PitNET is highly uneven. Graded from the
companion deep-research audit (`docs/deep-research-nf-pitnet-drug-lineage-validation.md`):

| Drug class | Representative agents | PitNET evidence | NFPA coverage | Evidence grade | Literature (PMID) |
|---|---|---|---|---|---|
| **HDAC** | panobinostat, vorinostat, CUDC-907 | Primary-cell HTS (9 PitNET); Cushing xenograft ↓65%, ACTH ↓77% | moderate–strong | **Deepest** | 38637883; 33000123 |
| **mTOR** | everolimus, AZD-8055 | GH3 + 6 human GH cultures; 4 prolactinomas 3/4 biochemical response | moderate (GH) | **Deepest** | 30066286; 37428396; 39736867 |
| **BET/BRD4** | JQ1, OTX015, ZBC-260 | GH3/MMQ G2/M arrest + apoptosis; **NFPA primary + xenograft** | **strong (NFPA models)** | Strong | 32246150; 31935194 |
| **MEK** | selumetinib, PD-0325901, binimetinib | corticotroph AtT20 + Cushing xenograft (TR4-dependent) | **corticotroph only** | Weak (lineage-limited) | 27708250 |
| **CDK4/6** | palbociclib | 1 NFPA regression (incidental, no lineage stratification) | **single case** | Lowest | 29852762 |
| EGFR/HER | lapatinib | phase-2a (NCT00939523): 0/4 met 40% shrinkage endpoint | failed example | Weak | 33150390 |

**Interpretation.** Three of five nominated classes (BET, HDAC, mTOR) have genuine PitNET
preclinical/clinical evidence, with HDAC(panobinostat) and mTOR(everolimus) the deepest and BET
the most complete in NFPA models — supporting the broad direction that the NFPT epigenetic/transcriptional
program is druggable. By contrast, **MEK evidence is confined to corticotroph/Cushing lineages and
CDK4/6 to a single NFPA case** — precisely the weakest-supported nominees — so these are demoted to
"low-confidence priorities pending NFPA-specific functional validation" in the main text. The
SF1-lineage methylation identity provides biological *plausibility* for epigenetic/kinase targeting
but is **not** itself evidence of drug sensitivity.

### S2.6 Limitations of the screen

1. Enrichr Combined Scores span orders of magnitude; absolute ranks are dominated by a few extreme
   BET instances (LINCS bromodomain bias). We therefore report at the **drug-class** level.
2. LINCS instances aggregate many cell lines, doses, and time points; the reversal score is a
   *signature-level* association, not NF-PitNET-cell-specific efficacy.
3. The framework is an Enrichr/LINCS approximation of classic CMap, not the `build02` rank-matrix.
4. Directionality was enforced by the four-query scheme, but wet-lab validation (organoids / primary
   cultures) is required before any therapeutic claim.
5. **No study has tested whether SF1-lineage status predicts response to any inhibitor**; the
   lineage→drug-sensitivity link is a hypothesis, not a prediction.

---

## References for Supplementary Methods

1. Shi L, et al. BRD4 as a therapeutic target for nonfunctioning and growth hormone pituitary adenoma. *Neuro Oncol*. 2020;22(8):1114-1125. PMID: 32246150.
2. Lines K, et al. Effects of epigenetic pathway inhibitors on corticotroph tumour AtT20 cells. *Endocr Relat Cancer*. 2020;27(3):163-174. PMID: 31935194.
3. Cheng Y, et al. Therapeutic potential of targeting Nrf2 by panobinostat in pituitary neuroendocrine tumors. *Acta Neuropathol Commun*. 2024. PMID: 38637883.
4. Zhang D, et al. Targeting Corticotroph HDAC and PI3-Kinase in Cushing Disease. *J Clin Endocrinol Metab*. 2021;106(1):e232-e246. PMID: 33000123.
5. Pivonello C, et al. Effect of combined treatment with a pan-PI3K inhibitor or an isoform-specific PI3K inhibitor and everolimus on cell proliferation in GH-secreting pituitary tumour. *Endocrine*. 2018. PMID: 30066286.
6. Lin AL, et al. The treatment of aggressive prolactinomas with everolimus. *Pituitary*. 2023. PMID: 37428396.
7. Yao Z, Chen H. Everolimus in pituitary tumor: a review of preclinical and clinical evidence. *Front Endocrinol (Lausanne)*. 2024. PMID: 39736867.
8. Zhang D, et al. Targeting the ERK pathway for the treatment of Cushing's disease. *Oncotarget*. 2016. PMID: 27708250.
9. Anderson E, et al. Regression of a nonfunctioning pituitary macroadenoma on the CDK4/6 inhibitor palbociclib: case report. *Neurosurg Focus*. 2018. PMID: 29852762.
10. Cooper O, et al. EGFR/ErbB2-Targeting Lapatinib Therapy for Aggressive Prolactinomas. *J Clin Endocrinol Metab*. 2021. PMID: 33150390.
