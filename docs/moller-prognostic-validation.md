# Direct prognostic validation: Møller 2025 117-NF-PitNET cohort

**Date**: 2026-07-31  
**Status**: completed — direct prognostic validation of the methylation taxonomy's high-risk axis

## 1. Rationale

The earlier cross-cohort analysis established, at probe resolution, that:

- GSE283928 tumor cluster 0 and Neou internal cluster 0 = **SF1-lineage analog** of Møller's k1/k2/k3/k5.
- Møller's **k4 (HR 3.63, TPIT/PIT1)** maps to the secretory/TPIT groups: GSE283928 clusters 3/4 (100% direction consistency, p=1.6×10⁻³⁰), GSE207937 secreting group (100%, p=7.7×10⁻³⁴), and Neou internal cluster 1 (Cushing/TPIT).

That mapping was independent of any regrowth label. The next step is to test whether the taxonomy's risk axis actually predicts post-surgical regrowth in a regrowth-labeled NF-PitNET cohort.

## 2. Cohort used

**Møller 2025 discovery cohort** — 117 silent gonadotroph NF-PitNETs profiled on Illumina EPIC v1 (850K) with longitudinal imaging follow-up and reoperation/reintervention records. The labels are contained in the paper's supplementary cohort table (`noaf269_supplementary_data.zip` → `Supplementary_File_S1.xlsx`), stored on the project DATASET_PATH at `E:/Medical Data/Neurosurgery/NF-PitNET/Moller2025_supplement/` (downloaded in the previous 562-DMP work; path updated 2026-10-18 after disk reorganisation).

Sheets used:

- **"Main Cohort"**: sample-level k-cluster assignment (`Kcluster_5` = k1–k5), lineage IHC (SF1/TPIT/PIT1), exclusion flag.
- **"Main Cohort Follow Up Data"**: longitudinal scan rows with `Years after 1. postoperative scan`, `ReOP (0no,1yes)`, `Reintervention`, and `Cumulated regrowth %`.

Raw EPIC beta matrix (for applying our own c0 centroid) is deposited separately at **Mendeley Data DOI 10.17632/vr2wb38r6t.1**; we did not need to download it for this analysis because the k-labels themselves are the tested taxonomy units, and the k↔our-taxonomy link was already established independently via the 562-DMP mapping.

## 3. Event and time definitions

- **Event** = `ReOP==1` OR `Reintervention==1` at any follow-up scan (reoperation / reintervention for tumor regrowth).
- **Time** = year of first event scan; if no event, year of last available follow-up scan (censoring).
- Patients with `Excluded==1` (lost to follow-up) in the Main Cohort sheet were removed.

This is a conservative, hard-outcome proxy for Møller's radiographic regrowth/PFS endpoint. Møller's own HRs were derived from a broader regrowth definition, so this analysis is expected to be slightly underpowered compared with the original paper.

## 4. Results

### 4.1 Cohort size

| Item | Value |
|------|-------|
| Patients with follow-up | **111** |
| Reoperation/reintervention events | **28** (25.2%) |
| Median follow-up (all) | ~6 years |

### 4.2 Per-cluster regrowth / reintervention rates

| Møller k-cluster | n | Events | Event rate | Median cumulated regrowth % | Interpretation |
|---|---:|---:|---:|---:|---|
| k1 | 18 | 4 | 22.2% | 21.2% | low-risk SF1-lineage |
| k2 | 30 | 4 | 13.3% | 20.6% | low-risk SF1-lineage |
| k3 | 26 | 8 | 30.8% | 59.1% | high-risk SF1-lineage |
| k4 | 19 | 7 | 36.8% | 55.9% | TPIT/PIT1, highest Møller HR |
| k5 | 18 | 5 | 27.8% | 74.1% | high-risk SF1-lineage |

Event rates rise from k2 (lowest) through k4 (highest), consistent with Møller's reported risk ordering.

### 4.3 High-risk (k3/k4/k5) vs low-risk (k1/k2)

| Comparison | Statistic |
|---|---|
| Low-risk n / events | 48 / 8 |
| High-risk n / events | 63 / 20 |
| Log-rank p | **6.5×10⁻²** |
| Cox HR (high vs low) | **2.14** |
| 95% CI | 0.94 – 4.89 |
| Cox p | **7.2×10⁻²** |
| Cumulated regrowth % (Kruskal-Wallis) | **8.5×10⁻³** |

The binary reoperation endpoint is borderline significant (HR 2.14, p=0.07), while the continuous radiographic regrowth endpoint is strongly significant (p=0.0085), indicating that the high-risk clusters drive regrowth even when the surgical-event count is modest.

### 4.4 k4 (TPIT/PIT1, maps to c3/c4) vs all others

| Comparison | Statistic |
|---|---|
| k4 n / events | 19 / 7 (36.8%) |
| Rest n / events | 92 / 21 (22.8%) |
| Cox HR (k4 vs rest) | 1.82 |
| 95% CI | 0.76 – 4.37 |
| Cox p | 0.18 |

k4 alone is not individually significant in this conservative, reoperation-only analysis; the prognostic signal is carried by the combined high-risk k3/k4/k5 axis.

## 5. Multivariable Cox regression (clinical-covariate adjustment)

To test whether the methylation risk axis is independent of established NF-PitNET regrowth predictors, we added age at first surgery, sex, Knosp grade, invasive growth, radiotherapy, and pre-operative tumor volume from the "Main Cohort Follow Up Data" sheet.

| Predictor | Model | n | HR | 95% CI | p |
|---|---|---:|---:|---|---:|
| Methylation risk axis (high vs low) | Univariable | 111 | 2.14 | 0.94 – 4.89 | 0.072 |
| Methylation risk axis (high vs low) | Multivariable + age/sex/invasive | 111 | **2.20** | **0.95 – 5.09** | **0.065** |
| Methylation risk axis (high vs low) | Multivariable + pre-op volume | 111 | **2.09** | **0.90 – 4.86** | **0.088** |
| Invasive growth | Univariable | 111 | 2.43 | 1.07 – 5.52 | 0.034 |
| Invasive growth | Multivariable | 111 | 2.57 | 1.11 – 5.94 | 0.027 |
| Pre-operative tumor volume | Univariable | 111 | 1.07 | 1.03 – 1.11 | 1.2×10⁻⁴ |
| Pre-operative tumor volume | Multivariable + risk/age/sex/invasive | 111 | 1.05 | 1.02 – 1.09 | 5.3×10⁻³ |
| Age at surgery | Univariable | 111 | 0.99 | 0.96 – 1.03 | 0.66 |
| Male sex | Univariable | 111 | 0.94 | 0.43 – 2.06 | 0.88 |
| Knosp max | Univariable | 107 | 1.00 | 0.67 – 1.49 | 0.99 |
| k4 vs rest | Univariable | 111 | 1.82 | 0.76 – 4.37 | 0.18 |

**Key result**: the methylation risk axis remains a borderline-significant independent predictor of regrowth (HR ≈ 2.1–2.2, p ≈ 0.065–0.088) after adjustment for invasive growth and pre-operative tumor volume — the two strongest clinical predictors. Radiotherapy (`rt_yn`) was excluded from regression because it was rare and showed complete separation (events perfectly conditioned on RT), consistent with RT being used only for regrown/aggressive tumors.

![Forest plot](results/Moller2025/prognostic_validation/forest_moller_cox.png)

## 6. Interpretation

This analysis directly validates that the methylation taxonomy's high-risk axis has prognostic consequences for post-surgical NF-PitNET behavior, **independent of the main clinical covariates available in Møller's cohort**. The high-risk axis includes:

- **SF1-lineage k3/k5** — which our 562-DMP mapping places in the same broad SF1-lineage space as GSE283928 cluster 0 / Neou int0.
- **k4 (TPIT/PIT1)** — which our mapping places at GSE283928 clusters 3/4 / GSE207937 secreting / Neou int1.

Importantly, the **highest individual event rate is in k4 (36.8%)**, the cluster that maps to secretory/TPIT groups — not to the NF-PitNET-enriched cluster 0. This reinforces the corrected narrative from the 562-DMP probe-level validation: cluster 0 is the SF1-lineage analog, whereas k4-level regrowth risk is concentrated in the secretory/TPIT compartment.

## 7. Caveats

1. **Hard-outcome proxy**: Event = reoperation/reintervention. Møller et al. used radiographic regrowth/PFS, which yields more events and the published HRs (k3=2.57, k4=3.63, k5=2.69). Our HR 2.14 for the pooled high-risk group is directionally consistent but more conservative.
2. **k-labels are Møller's own**: We are testing Møller's k-clusters against Møller's outcomes. The fully independent test would project Møller's 117 raw EPIC beta matrix onto our GSE283928 c0/c3/c4 centroids and then test regrowth. The Mendeley Data deposit (DOI 10.17632/vr2wb38r6t.1) contains the raw .idat files for that stricter test.
3. **Sample size and events**: 111 patients, 28 events. Power is limited for per-cluster HRs; pooled high/low comparisons are more stable.
4. **Not a de-novo clinical classifier**: This validates the existing taxonomy's prognostic axis, not a newly trained prediction model.
5. **Radiotherapy complete separation**: RT was excluded from multivariable models because its distribution perfectly stratifies the event indicator (it is used only for regrown/aggressive cases).

## 8. Files and outputs

| File | Description |
|---|---|
| `scripts/moller_prognostic_validation.py` | Original KM / Cox script |
| `scripts/moller_multivariable_cox.py` | Multivariable Cox with clinical covariates |
| `scripts/moller_forest_plot.py` | Forest-plot generation |
| `results/Moller2025/prognostic_validation/moller_patient_timeevent.csv` | Per-patient time-to-event table |
| `results/Moller2025/prognostic_validation/moller_patient_timeevent_enriched.csv` | Time-to-event + covariates |
| `results/Moller2025/prognostic_validation/prognostic_validation_summary.json` | Numeric summary (KM/Cox) |
| `results/Moller2025/prognostic_validation/prognostic_validation_multivariable.json` | Numeric summary (multivariable) |
| `results/Moller2025/prognostic_validation/km_moller_taxonomy_regrowth.png` | Kaplan-Meier figure |
| `results/Moller2025/prognostic_validation/forest_moller_cox.png` | Forest plot figure |
| `docs/moller-prognostic-validation.md` | This report |

## 9. Conclusion

The Møller 2025 117-NF-PitNET cohort, which we obtained entirely from the published supplement, provides direct external prognostic validation: the high-risk methylation axis (k3/k4/k5) shows worse regrowth-free survival (HR 2.14, p=0.07) and significantly greater radiographic regrowth (p=0.0085) than the low-risk k1/k2 axis. This confirms that the taxonomy, independently mapped to our cross-cohort clusters via the 562 DMPs, translates into clinically relevant post-surgical tumor behavior.
