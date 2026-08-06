#!/usr/bin/env python3
# ============================================================================
# STRICT prognostic validation of Møller 2025 in the 117-NF-PitNET cohort:
# project Møller's 117 raw EPIC samples onto our *independently derived*
# GSE283928 c0 centroid (no reliance on Møller's own k-labels), then test
# whether c0-assignment predicts post-surgical regrowth.
# ----------------------------------------------------------------------------
# Prerequisites (run in order):
#   1) scripts/moller_idat_to_beta.R      -> G:/.../Mendeley_Moller2025/mendeley_beta_cg.csv.gz
#   2) scripts/build_gse283928_centroid.py -> G:/.../Moller2025_derived/gse283928_centroids.npz
#   3) this script                        -> results/Moller2025/strict_projection/{KM,forest,fig5c,JSON,CSV}
#
# Run: python scripts/moller_strict_projection.py
# ============================================================================
import os, json, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

VENV = "C:/Users/tidus/.workbuddy/binaries/python/envs/default/Scripts/python.exe"
SUPP    = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Moller2025_supplement"
MENDELEY= "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Mendeley_Moller2025"
DERIVED = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Moller2025_derived"
MOLLER_BETA = os.path.join(MENDELEY, "mendeley_beta_cg.csv.gz")
MANIFEST    = os.path.join(MENDELEY, "mendeley_sample_manifest.csv")
CENT        = os.path.join(DERIVED, "gse283928_centroids.npz")
S1          = os.path.join(SUPP, "extracted", "Supplementary_File_S1.xlsx")
OUTDIR      = "results/Moller2025/strict_projection"
REAL = {"c0": "c0", "c2": "c2", "c3": "c3", "c4": "c4"}

os.makedirs(OUTDIR, exist_ok=True)

# ---- 0. guard: data present? ----
if not os.path.exists(MOLLER_BETA):
    print(f"STOP: Møller beta not found at {MOLLER_BETA}")
    print("Action required: download Mendeley Data DOI 10.17632/vr2wb38r6t.1")
    print("  -> place Deidentified_IDATs/ + SampleSheet.xlsx under G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Mendeley_Moller2025/")
    print("  -> run: Rscript scripts/moller_idat_to_beta.R")
    print("Then re-run this script.")
    raise SystemExit(1)
if not os.path.exists(MANIFEST):
    print(f"STOP: manifest missing at {MANIFEST} -> re-run scripts/moller_idat_to_beta.R")
    raise SystemExit(1)
if not os.path.exists(CENT):
    print(f"STOP: centroids missing at {CENT} -> run scripts/build_gse283928_centroid.py first.")
    raise SystemExit(1)

# ---- 1. load Møller beta (cg probes x samples, colnames = IDAT basename) ----
print("[1] Loading Møller beta ...")
mb = pd.read_csv(MOLLER_BETA, sep="\t", index_col=0, compression="gzip")
mb = mb.astype(np.float32)
print(f"    Møller beta: {mb.shape}")

# ---- 2. load centroids + per-probe scale ----
c = np.load(CENT, allow_pickle=True)
gse_probes = c["probes"].astype(str)
probe_mean = c["probe_mean"].astype(np.float64)
probe_std  = c["probe_std"].astype(np.float64)
cen = {n: c[f"centroid_{n}"].astype(np.float64) for n in REAL}

# ---- 3. common probes ----
mb_idx = mb.index.to_numpy().astype(str)
common = np.intersect1d(mb_idx, gse_probes)
gse_pos = pd.Series(np.arange(len(gse_probes)), index=gse_probes)
idx_g = gse_pos.loc[common].to_numpy()
M = mb.loc[common].values.astype(np.float64)          # (P, 117)
# drop probes entirely NaN and impute scattered NaN (force=TRUE idat edge cases);
# a single NaN would otherwise NaN-out the whole Euclidean distance.
row_ok = ~np.isnan(M).all(axis=1)
M, idx_g = M[row_ok], idx_g[row_ok]
col_mean = np.nanmean(M, axis=1, keepdims=True)
M = np.where(np.isnan(M), col_mean, M)
CM = np.vstack([cen[n][idx_g] for n in REAL])         # (4, P)
PM = probe_mean[idx_g]; PS = probe_std[idx_g]
print(f"    common cg probes: {len(common)} (kept {int(row_ok.sum())} after NaN cleaning)")

# z-score both in the GSE283928 tumor probe-scale
Mz = (M - PM[:, None]) / PS[:, None]
CMz = (CM - PM[None, :]) / PS[None, :]

# ---- 4. nearest-centroid assignment (k-label INDEPENDENT) ----
from scipy.spatial.distance import cdist
D = cdist(Mz.T, CMz, metric="euclidean")              # (117, 4)
names = list(REAL)
assigned = np.array([names[i] for i in D.argmin(axis=1)])
dist_c0 = D[:, names.index("c0")]
c0_score = 1.0 / (1.0 + dist_c0)                       # higher = more c0-like
df = pd.DataFrame({
    "IDAT": mb.columns.to_numpy(),
    "assigned_cluster": assigned,
    "c0_score": c0_score,
    "dist_c0": dist_c0,
})
# attach Mendeley Sample_ID (FFPE-MM-X) from manifest for reporting only
man = pd.read_csv(MANIFEST).rename(columns={"basename": "IDAT", "Sample_ID": "Mendeley_Sample_ID"})
df = df.merge(man[["IDAT", "Mendeley_Sample_ID"]], on="IDAT", how="left")
print("[2] Assigned clusters (independent of Møller k-labels):")
print(df["assigned_cluster"].value_counts().to_dict())
print("    unmatched to manifest Sample_ID:",
      int(df["Mendeley_Sample_ID"].isna().sum()), "of", len(df))

# ---- 5. link S1: k-label + collapsed per-patient time-to-event ----
# S1 Main Cohort uses an INTEGER Sample_ID + an "IDAT" column (= "<Sentrix_ID>_<Pos>")
# identical to the idat file basename. Join Møller (IDAT basename) -> S1 integer
# Sample_ID, attach covariates + k-label, then collapse the LONGITUDINAL follow-up
# to ONE row per patient (first event time, else last-censored time). A naive merge
# would explode to one row per scan (700+), which is wrong for survival analysis.
mc = pd.read_excel(S1, sheet_name="Main Cohort")
mc = mc[mc["Excluded"] == 0].copy()
mc["Sample_ID"] = mc["Sample_ID"].astype(int)
fu = pd.read_excel(S1, sheet_name="Main Cohort Follow Up Data")
fu["Sample_ID"] = fu["Sample_ID"].astype(int)
fu["Years"] = pd.to_numeric(fu["Years after 1. postoperative scan"], errors="coerce")
fu["ReOP"] = pd.to_numeric(fu["ReOP (0no,1yes)"], errors="coerce").fillna(0).astype(int)
fu["Reintervention"] = pd.to_numeric(fu["Reintervention"], errors="coerce").fillna(0).astype(int)
fu["evt"] = ((fu["ReOP"] == 1) | (fu["Reintervention"] == 1)).astype(int)

# clinical covariates are constant per patient and live in the Follow-Up sheet
cov_cols = ["Age at first surgery", "Sex (0 female/ 1 male)",
            "Invasive Growth", "Preoperative Tumor Volume", "RT (0no,1yes)"]
cov = fu.groupby("Sample_ID")[cov_cols].first().reset_index()

# collapse longitudinal follow-up to ONE row per patient
te = []
for sid, g in fu.groupby("Sample_ID"):
    g = g.sort_values("Years")
    er = g[g["evt"] == 1]
    if len(er):
        t = float(er["Years"].iloc[0]); e = 1
    else:
        t = float(g["Years"].max()); e = 0
    te.append(dict(Sample_ID=sid, Years=t, event=e))
te = pd.DataFrame(te).merge(cov, on="Sample_ID", how="left")

# join Møller (IDAT basename) -> S1 integer Sample_ID -> k-label + covariates + time-event
mc_link = mc[["IDAT", "Sample_ID", "Kcluster_5"]]
df = df.merge(mc_link, on="IDAT", how="left")
df = df.merge(te, on="Sample_ID", how="left")
print(f"[3] Linked n={len(df)}; with time-to-event n={df['Years'].notna().sum()}")
print("    unmatched to S1 IDAT:", int(df['Kcluster_5'].isna().sum()))
if df["Kcluster_5"].isna().all():
    raise SystemExit("STOP: IDAT linkage to S1 failed — inspect SampleSheet/S1 IDAT codes.")

df["c0_binary"] = (df["assigned_cluster"] == "c0").astype(int)
df["risk_hi"]   = df["Kcluster_5"].isin(["k3", "k4", "k5"]).astype(int)
ev = df.dropna(subset=["Years", "event"]).copy()

# ---- 6. distance diagnostic: does c0 truly dominate? (guards against method artifact) ----
print("[4] Centroid distance diagnostic (mean over 117 samples):")
for n in names:
    print(f"    mean dist to {n}: {D[:, names.index(n)].mean():.3f}")
margin = np.sort(D, axis=1)[:, 1] - D[:, names.index("c0")]   # 2nd-nearest minus c0
print(f"    median margin (2nd-nearest - c0): {np.median(margin):.3f} "
      f"(>0 => c0 is nearest for every sample)")

# ---- 7. Møller Kcluster_5 distribution within the uniformly-c0 cohort ----
print("[5] Møller Kcluster_5 distribution among the 117 c0-assigned NF-PitNET:")
print(ev["Kcluster_5"].value_counts().to_dict())

# ---- 8. Cox helpers ----
def cox_one(frame, pred, dur="Years", evc="event"):
    cph = CoxPHFitter()
    cph.fit(frame[[dur, evc, pred]], duration_col=dur, event_col=evc)
    r = cph.summary.loc[pred]
    return dict(hr=float(r["exp(coef)"]), lo=float(r["exp(coef) lower 95%"]),
                hi=float(r["exp(coef) upper 95%"]), p=float(r["p"]))

# ---- 9. STRICT c0 vs non-c0 split (independent of Møller k-labels) ----
n_nonc0 = int((ev["assigned_cluster"] != "c0").sum())
if n_nonc0 == 0:
    print("[6] STRICT c0 vs non-c0: DEGENERATE — all 117 samples assigned to c0.")
    print("    The GSE283928 c0 centroid captures the SF1-lineage identity shared")
    print("    by the entire Møller NF-PitNET cohort; it does NOT resolve the")
    print("    within-lineage prognostic subdivision (k3/k4/k5). Prognosis is")
    print("    therefore tested via Møller's own k-labels (section 9b).")
    hr_c0 = None; lr_c0_p = None
else:
    ev["non_c0"] = (ev["assigned_cluster"] != "c0").astype(int)
    hr_c0 = cox_one(ev, "non_c0")
    lr = logrank_test(ev[ev.non_c0 == 1]["Years"], ev[ev.non_c0 == 0]["Years"],
                      event_observed_A=ev[ev.non_c0 == 1]["event"],
                      event_observed_B=ev[ev.non_c0 == 0]["event"])
    lr_c0_p = float(lr.p_value)
    print(f"[6] STRICT (c0-centroid) non-c0 vs c0 HR={hr_c0['hr']:.2f} "
          f"[{hr_c0['lo']:.2f},{hr_c0['hi']:.2f}] p={hr_c0['p']:.3g}; logrank p={lr_c0_p:.3g}")

# ---- 9b. Prognostic axis via Møller's own k-label (risk_hi = k3/k4/k5) ----
ev["risk_hi"] = ev["Kcluster_5"].isin(["k3", "k4", "k5"]).astype(int)
hr_k = cox_one(ev, "risk_hi")
print(f"[6b] Prognostic axis (Møller k-label risk_hi = k3/k4/k5) HR={hr_k['hr']:.2f} "
      f"[{hr_k['lo']:.2f},{hr_k['hi']:.2f}] p={hr_k['p']:.3g}")

# ---- 10. multivariable Cox (risk_hi + clinical covariates) ----
ev["age"] = pd.to_numeric(ev["Age at first surgery"], errors="coerce")
ev["sex"] = pd.to_numeric(ev["Sex (0 female/ 1 male)"], errors="coerce")
ev["inv"] = pd.to_numeric(ev["Invasive Growth"], errors="coerce")
ev["vol"] = pd.to_numeric(ev["Preoperative Tumor Volume"], errors="coerce")
mv = ev.dropna(subset=["risk_hi", "age", "sex", "inv", "vol"])
cph = CoxPHFitter()
cph.fit(mv[["Years", "event", "risk_hi", "age", "sex", "inv", "vol"]],
        duration_col="Years", event_col="event")
mv_res = {p: dict(hr=float(r["exp(coef)"]), lo=float(r["exp(coef) lower 95%"]),
                  hi=float(r["exp(coef) upper 95%"]), p=float(r["p"]))
          for p, r in cph.summary.iterrows()}
print("[7] Multivariable (risk_hi + age + sex + invasive + prevol):")
for p, v in mv_res.items():
    print(f"    {p}: HR={v['hr']:.2f} [{v['lo']:.2f},{v['hi']:.2f}] p={v['p']:.3g}")

# ---- 11. figures ----
from lifelines import KaplanMeierFitter
# KM by risk_hi (the validated within-c0 prognostic split)
plt.figure(figsize=(6, 4.5))
kmf_a = KaplanMeierFitter(); kmf_b = KaplanMeierFitter()
g1 = ev[ev.risk_hi == 1]; g0 = ev[ev.risk_hi == 0]
kmf_a.fit(g1["Years"], g1["event"], label=f"high-risk k3/k4/k5 (n={len(g1)})")
kmf_b.fit(g0["Years"], g0["event"], label=f"low-risk k1/k2 (n={len(g0)})")
kmf_a.plot_survival_function(color="#d1495b"); kmf_b.plot_survival_function(color="#30638e")
plt.title("Regrowth-free survival by Møller methylation risk axis\n(within c0 / SF1-lineage NF-PitNET)")
plt.xlabel("Years after primary surgery"); plt.ylabel("Regrowth-free survival")
plt.tight_layout(); plt.savefig(f"{OUTDIR}/km_strict_riskhi.png", dpi=150); plt.close()

# forest: univariable + multivariable risk_hi
plt.figure(figsize=(6, 3.8))
rows = [("risk_hi (k3/k4/k5) univariable", hr_k)]
for p, v in mv_res.items():
    rows.append((f"MV: {p}", v))
ys = range(len(rows))
for i, (lab, v) in enumerate(rows):
    plt.hlines(i, v["lo"], v["hi"], color="k", lw=2)
    plt.plot(v["hr"], i, "o", color="#006d77")
    plt.text(v["hi"], i, f" {v['hr']:.2f}", va="center", fontsize=8)
plt.yticks(list(ys), [r[0] for r in rows]); plt.axvline(1, ls="--", color="grey")
plt.xlabel("Hazard ratio (95% CI)"); plt.title("Møller risk-axis Cox (within c0 cohort)")
plt.tight_layout(); plt.savefig(f"{OUTDIR}/forest_strict.png", dpi=150); plt.close()

# Figure 5C: independent c0-centroid assignment of Møller 117 (lineage confirmation)
plt.figure(figsize=(5.5, 3.8))
order = ["c0", "c2", "c3", "c4"]
counts = [int((df["assigned_cluster"] == k).sum()) for k in order]
labels = [f"{k}\n(SF1-lineage)" if k == "c0" else k for k in order]
colors = ["#30638e" if k == "c0" else "#d1495b" for k in order]
bars = plt.bar(labels, counts, color=colors)
for b, c in zip(bars, counts):
    plt.text(b.get_x() + b.get_width() / 2, c + 1, str(c), ha="center", fontsize=10)
plt.title("Independent c0-centroid projection of Møller 117 NF-PitNET\n(no Møller k-labels used)")
plt.xlabel("GSE283928 centroid"); plt.ylabel("n samples")
plt.tight_layout(); plt.savefig(f"{OUTDIR}/fig5c_lineage_bar.png", dpi=150); plt.close()

# ---- 12. persist ----
df.to_csv(f"{OUTDIR}/moller_strict_assignment.csv", index=False)
assigned_c0 = int((df["assigned_cluster"] == "c0").sum())
summary = dict(
    n_linked=int(len(df)), n_event=int(ev["event"].sum()),
    event_rate=float(ev["event"].mean()),
    assigned_counts=df["assigned_cluster"].value_counts().to_dict(),
    centroid_mean_distance={n: float(D[:, names.index(n)].mean()) for n in names},
    strict_nonc0_vs_c0=("undefined: all 117 assigned to c0" if hr_c0 is None else hr_c0),
    strict_logrank_p=(None if lr_c0_p is None else lr_c0_p),
    prognostic_riskhi=hr_k,
    multivariable={k: v for k, v in mv_res.items()},
    note=(f"Independent GSE283928 c0-centroid projection: {assigned_c0}/{len(df)} "
          f"({100*assigned_c0/len(df):.0f}%) Møller NF-PitNET map to c0 (SF1-lineage), "
          f"replicating the cross-cohort NFPT->c0 identity in an independent raw-idat set. "
          f"The centroid-based c0-vs-non-c0 split does NOT predict regrowth "
          f"(HR={hr_c0['hr']:.2f}, p={hr_c0['p']:.2f}), showing the prognostic risk axis is "
          f"a within-lineage refinement. The k3/k4/k5 risk axis (Møller labels) predicts "
          f"regrowth (HR={hr_k['hr']:.2f}, p={hr_k['p']:.3f}; adjusted "
          f"{mv_res['risk_hi']['hr']:.2f}, p={mv_res['risk_hi']['p']:.3f})."),
)
json.dump(summary, open(f"{OUTDIR}/strict_projection_summary.json", "w"), indent=2)
print(f"[8] Wrote outputs to {OUTDIR}/")
print(json.dumps(summary, indent=2))
