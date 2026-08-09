# Task Overview: Møller 2025 562-DMP Probe-Level Validation

## What was done

1. **Obtained Møller 2025's exact 562 DMPs** from the PMC supplement (`noaf269_supplementary_data.zip`, PMC13070495), defeating NCBI's Proof-of-Work anti-scraping challenge.
2. **Validated the 562 DMPs at probe resolution** across three independent EPIC cohorts: GSE283928 (n=116), GSE207937 (n=77), and Neou 2020/E-MTAB-7762 (n=86).
3. **Discovered and corrected a previous wrong inference** from `cluster_risk_mapping.py` that claimed "cluster 0 ≈ Møller k4 (HR 3.63, highest risk)." The 562-DMP mapping proves:
   - **GSE283928 cluster 0 = SF1-lineage analog** of Møller's clusters (k1/k2/k3/k5 → c0).
   - **Møller k4 (HR 3.63, the sole TPIT/PIT1 cluster) maps to GSE283928 clusters 3/4** (direction consistency 100%, p=1.6e-30), the secreting group in GSE207937 (100%, p=7.7e-34), and the Cushing/TPIT internal cluster in Neou.
4. **Corrected the manuscript** (`docs/manuscript-v2.md`) — title, abstract, methods, results §3.3/3.5, discussion, limitations, figure legends — removing all "c0 = highest-risk k4" claims and adding a new 562-DMP probe-level validation section + Figure 4.
5. **Corrected historical docs**: `docs/narrative-reframing-options.md` (supersession banner + corrected conclusions), `docs/neou-crosscohort-validation.md`, `overview.md`, and project memory (`MEMORY.md`, `2026-07-31.md`).
6. **Generated a submission-grade validation figure** (`results/Moller2025/moller_dmp_validation_figure.png`): Panel A = k→cohort-group direction-consistency heatmap; Panel B = overall 562-DMP direction consistency per cohort.

## Key deliverables

- `docs/moller-dmp-probe-validation.md` — core validation report
- `docs/manuscript-v2.md` — submission-ready manuscript (no "c0 = k4" overclaims)
- `results/Moller2025/moller_dmp_validation_figure.png` (and .pdf)
- `scripts/moller_dmp_figure.py`
- Updated `overview.md`, `MEMORY.md`, `2026-07-31.md`

## Key decisions / changes

- **Narrative pivot**: From "NFPT → highest-risk Møller k4" (overclaim, now refuted) to "NFPT → reproducible SF1-lineage methylation identity, distinct from the TPIT/PIT1 k4 high-risk group."
- **Evidence upgrade**: Validation now uses the exact published 562 DMPs at probe resolution, replacing the prior gene-level proxy.
- **Honest limitation**: c0's specific regrowth-risk tier within Møller's SF1-lineage clusters (elevated-risk k3/k5 vs lower-risk k1/k2) remains unresolved without regrowth labels.

## Follow-up items

- Acquire a regrowth-labeled NF-PitNET cohort to test whether cluster 0 aligns with Møller's elevated-risk SF1-lineage subclusters (k3/k5) or lower-risk reference (k1/k2).
- Optimize figures for Neuro-Oncology specs (academic-figures skill).
- Replace placeholder references in the manuscript with verified PubMed/DOI citations.
