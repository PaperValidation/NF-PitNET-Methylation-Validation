# NF-PitNET Methylation Taxonomy: Cross-Cohort Validation

**SF1-lineage identity, not the highest-risk k4 cluster, defines non-functioning pituitary neuroendocrine tumors: independent cross-cohort and cross-platform methylation validation**

This repository contains the full analysis code and processed results for an independent, cross-cohort validation of the NF-PitNET methylation taxonomy (target journal: *Neuro-Oncology*).

## Overview

NF-PitNETs reproducibly map to the GSE283928 tumor cluster 0 (c0, an SF1-lineage methylation identity) across four independent cohorts spanning two microarray platforms (mean 84.3% c0 mapping; Neou 81.2%, GSE207937 75.6%, Moller-2025 raw-idat 86.3%, GSE115783 HM450 94.1%). Moller's highest-risk k4 cluster (TPIT/PIT1, HR 3.63) maps to the secretory clusters 3/4, not to c0 (probe-level p = 1.6×10⁻³⁰), correcting the common conflation of "highest-risk cluster" with "aggressive NF-PitNET." The within-lineage k-subdivision is recoverable from independent minfi processing (leave-one-out 99.1% nearest-centroid / 93.7% random forest), and the k3/k4/k5 high-risk axis shows a trend toward increased regrowth risk (HR 2.14, p = 0.072; adjusted HR 2.09, p = 0.088).

## Data sources

All raw data are third-party and publicly available; **no raw data are redistributed here** (only processed results and code). See [Data Availability](docs/manuscript-v1.md) for the full statement:

| Cohort | Platform | N (composition) | Source | Accession |
|---|---|---|---|---|
| Neou 2020 | EPIC-850K | 86 (32 NF-PitNET; 54 functioning) | ArrayExpress/BioStudies | E-MTAB-7762 |
| GSE283928 | EPIC-850K | 118 (4 tumor clusters: c0 n=70, c2 n=2, c3 n=22, c4 n=24) | GEO (Belakhoua et al., 2025) | GSE283928 |
| GSE207937 | EPIC-850K | 77 (45 NS/NF-PitNET; 32 functioning) | GEO (Silva-Junior et al., 2023) | GSE207937 |
| GSE115783 | HM450 | 40 (34 NFPA tumors; 6 healthy-pituitary controls) | GEO (Kober/Bujko et al., 2018/2019) | GSE115783 |
| Moller 2025 idat | EPIC-850K | 117 (silent gonadotroph NF-PitNET, k1-k5, longitudinal regrowth follow-up) | Mendeley Data | DOI 10.17632/vr2wb38r6t.1 |
| GSE209903 | RNA (transcriptome) | differential expression (NF-PitNET vs functioning) | GEO (Silva-Junior et al., 2023) | GSE209903 |

**Roles:** GSE283928 supplies the reference tumor clusters (c0/c2/c3/c4) used as nearest-centroid targets; Neou and GSE207937 are independent EPIC-850K validation cohorts; GSE115783 is the cross-platform (HM450) external validation cohort with embedded healthy-pituitary controls; the Moller-2025 raw idat are used for the fully de novo, k-label-free lineage projection and the direct prognostic test; GSE209903 provides the transcriptomic concordance check.

## Repository layout

```
scripts/    All Python (52) and R (16) analysis scripts, fully runnable
results/    Processed beta matrices, centroids, permutation outputs, figures, tables
docs/       Manuscript draft (manuscript-v1.md), supplementary materials, review reports
experience-archive/  Documented analysis pitfalls and reusable workflows
```

## Reproducing the analysis

1. Download raw data from the sources above into a local `DATASET_PATH` (paths are read from a local `.env`, which is **not** committed).
2. Follow `docs/supplementary-methods.md` for software versions and random seeds (minfi v1.48.0 / R 4.3.3; scikit-learn 1.5.1, random_state=42; 20,000 permutations, seed=0).
3. Run scripts in dependency order (see headers of each script; key pipeline entry points: `moller_idat_to_beta.R`, `neou_crosscohort.py`, `nn_mixing_permutation.py`, `build_supplementary_tables.py`).

## Key results files

- `results/Neou/joint_pca_coords.csv` — 281-sample joint PCA (PC1-20)
- `results/Neou/nn_mixing_permutation.json` — NN-mixing label-permutation p-value (p < 0.001)
- `results/Moller2025/strict_projection/` — 117-idat de novo projection (86.3% → c0)
- `results/Moller2025/within_lineage_562dmp/` — within-lineage k recovery (99.1%/93.7% LOO)
- `results/Moller2025/prognostic_validation/` — Cox regrowth models (HR 2.14 / adjusted 2.09)
- `results/Neou/supplementary_table1-3_*.csv` — supplementary tables
- `results/Moller2025/figure4_merged.png` — merged multi-panel figure

## AI disclosure

AI-assisted tools (WorkBuddy with large language models) were used as productivity aids for: (i) writing and debugging the analysis scripts in `scripts/`; (ii) drafting and copyediting the manuscript and supplementary materials; and (iii) literature and data-source verification. All scientific content, data interpretation, statistical analyses, and conclusions were generated, reviewed, and validated by the authors. The AI tools did not make independent scientific judgments, and all code and results are publicly available for independent verification. This statement mirrors the manuscript's AI Disclosure Statement (see `docs/manuscript-v1.md`).

## Ethics

This is a secondary analysis of publicly available, de-identified genomic datasets obtained from GEO, ArrayExpress/BioStudies, and Mendeley Data. No new human-subjects data were generated; the original studies obtained ethical approval and informed consent as described in their respective publications. No additional ethics approval was required for this in-silico reanalysis.

## License

Analysis code: for academic use; see manuscript for citation. Raw data remain under the original sources' terms.

---
*This repository is under active development. Versioned DOI will be minted via Zenodo upon acceptance.*
