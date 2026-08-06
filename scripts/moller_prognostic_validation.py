#!/usr/bin/env python
"""
Direct prognostic validation of the NF-PitNET methylation taxonomy in Møller 2025's
independent 117-NF-PitNET cohort (EPIC 850K, longitudinal imaging follow-up).

Data sources (already on disk, from noaf269_supplementary_data.zip):
  - Supplementary_File_S1.xlsx :: "Main Cohort"           -> Sample_ID, Kcluster_5 (k1-k5), Excluded
  - Supplementary_File_S1.xlsx :: "Main Cohort Follow Up Data" -> longitudinal regrowth / ReOP / Reintervention

We test whether the methylation taxonomy (k-clusters, which we independently mapped to our
c0/SF1-lineage vs k4/TPIT-PIT1 framework via the 562 DMPs) predicts post-surgical
regrowth / reintervention in an external NF-PitNET cohort.

Event definition: ReOP (reoperation) OR Reintervention == 1 at any follow-up scan.
Time: years after first post-operative scan (event time if event, else last-censored time).
"""
import openpyxl, json, os, numpy as np, pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats

ROOT = "results/Moller2025"   # derived outputs stay alongside code on C:
SUPP = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Moller2025_supplement"
S1 = os.path.join(SUPP, "extracted", "Supplementary_File_S1.xlsx")
OUT = os.path.join(ROOT, "prognostic_validation")
os.makedirs(OUT, exist_ok=True)

# ---- load Main Cohort (k-labels) ----
mc = pd.read_excel(S1, sheet_name="Main Cohort")
mc = mc[mc["Excluded"] == 0].copy()            # drop lost-to-follow-up
mc["Sample_ID"] = mc["Sample_ID"].astype(int)
mc = mc[["Sample_ID", "Kcluster_5", "SF1", "Tpit", "PIT1", "Lineage", "Gender"]]

# ---- load Follow Up Data (longitudinal regrowth) ----
fu = pd.read_excel(S1, sheet_name="Main Cohort Follow Up Data")
fu["Sample_ID"] = fu["Sample_ID"].astype(int)
fu["Years"] = pd.to_numeric(fu["Years after 1. postoperative scan"], errors="coerce")
fu["ReOP"] = pd.to_numeric(fu["ReOP (0no,1yes)"], errors="coerce").fillna(0).astype(int)
fu["Reintervention"] = pd.to_numeric(fu["Reintervention"], errors="coerce").fillna(0).astype(int)
fu["RegrowthPct"] = pd.to_numeric(fu["Cumulated regrowth %"], errors="coerce")
fu["Event"] = ((fu["ReOP"] == 1) | (fu["Reintervention"] == 1)).astype(int)

# ---- collapse to per-patient time-to-event ----
rows = []
for sid, g in fu.groupby("Sample_ID"):
    g = g.sort_values("Years")
    ev_rows = g[g["Event"] == 1]
    if len(ev_rows):
        t = float(ev_rows["Years"].iloc[0]); ev = 1
    else:
        t = float(g["Years"].max()); ev = 0
    final_regrow = float(g["RegrowthPct"].max()) if g["RegrowthPct"].notna().any() else np.nan
    rows.append(dict(Sample_ID=sid, time=t, event=ev, final_regrow_pct=final_regrow))
pt = pd.DataFrame(rows)
df = pt.merge(mc, on="Sample_ID", how="inner")
print(f"Patients with follow-up: {len(df)}  | events (ReOP/Reintervention): {df['event'].sum()} "
      f"({100*df['event'].mean():.1f}%)")

# ---- taxonomy grouping (our independent 562-DMP mapping) ----
# k1,k2,k3,k5 = SF1-lineage  -> our c0 axis
# k4          = TPIT/PIT1    -> our c3/c4 axis (Møller highest-risk HR 3.63)
df["group"] = df["Kcluster_5"].map(
    lambda k: "k4_TPITPIT1" if k == "k4" else "SF1lineage_k1k2k3k5")
# also Møller's own high/low split
df["risk"] = df["Kcluster_5"].map(
    lambda k: "high(k3k4k5)" if k in ("k3", "k4", "k5") else "low(k1k2)")

# ---- per-k descriptive ----
print("\n=== Per-kcluster regrowth/reintervention ===")
tbl = (df.groupby("Kcluster_5")
         .agg(n=("event", "size"), events=("event", "sum"),
              median_fu=("time", "median"),
              median_regrow=("final_regrow_pct", "median"))
         .reindex(["k1", "k2", "k3", "k4", "k5"]))
tbl["event_rate_%"] = (100 * tbl["events"] / tbl["n"]).round(1)
print(tbl.to_string())

# ---- (a) high vs low risk: log-rank + Cox ----
def cox_hr(frame, predictor_col, reference_note=""):
    cph = CoxPHFitter()
    cph.fit(frame[["time", "event", predictor_col]], duration_col="time", event_col="event")
    row = cph.summary.loc[predictor_col]
    return dict(hr=float(row["exp(coef)"]),
                lo=float(row["exp(coef) lower 95%"]),
                hi=float(row["exp(coef) upper 95%"]),
                p=float(row["p"]))

hi = df[df["risk"] == "high(k3k4k5)"]
lo = df[df["risk"] == "low(k1k2)"]
lr = logrank_test(hi["time"], lo["time"], event_observed_A=hi["event"], event_observed_B=lo["event"])
print(f"\n=== High(k3k4k5) vs Low(k1k2) ===")
print(f"Log-rank p = {lr.p_value:.3e}")
hr = cox_hr(df.assign(_x=(df["risk"] == "high(k3k4k5)").astype(int)), "_x")
print(f"Cox HR (high vs low) = {hr['hr']:.2f} 95% CI [{hr['lo']:.2f}, {hr['hi']:.2f}]  p = {hr['p']:.3e}")

# ---- (b) k4 vs rest (our c3/c4 axis) ----
k4 = df[df["Kcluster_5"] == "k4"]; rest = df[df["Kcluster_5"] != "k4"]
lr4 = logrank_test(k4["time"], rest["time"], event_observed_A=k4["event"], event_observed_B=rest["event"])
print(f"\n=== k4(TPIT/PIT1) vs rest ===  Log-rank p = {lr4.p_value:.3e}")
hr4 = cox_hr(df.assign(_x=(df["Kcluster_5"] == "k4").astype(int)), "_x")
print(f"Cox HR (k4 vs rest) = {hr4['hr']:.2f} 95% CI [{hr4['lo']:.2f}, {hr4['hi']:.2f}]  p = {hr4['p']:.3e}")

# ---- (c) continuous regrowth % by high/low ----
g_hi = df.loc[df["risk"] == "high(k3k4k5)", "final_regrow_pct"].dropna()
g_lo = df.loc[df["risk"] == "low(k1k2)", "final_regrow_pct"].dropna()
kw = stats.kruskal(g_hi, g_lo)
print(f"\n=== Cumulated regrowth %: high vs low Kruskal-Wallis p = {kw.pvalue:.3e} ===")

# ---- KM curves (figure) ----
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 11, "figure.dpi": 200})
fig, ax = plt.subplots(figsize=(6.2, 4.6))
km_lo = KaplanMeierFitter(); km_hi = KaplanMeierFitter()
km_lo.fit(lo["time"], lo["event"], label=f"Low-risk k1/k2 (n={len(lo)}, ev={lo['event'].sum()})")
km_hi.fit(hi["time"], hi["event"], label=f"High-risk k3/k4/k5 (n={len(hi)}, ev={hi['event'].sum()})")
km_lo.plot_survival_function(ax=ax, color="#2c7fb8", lw=2)
km_hi.plot_survival_function(ax=ax, color="#d7301f", lw=2)
ax.set_xlabel("Years after primary surgery"); ax.set_ylabel("Regrowth-free survival")
ax.set_title("Møller 117 NF-PitNET: methylation taxonomy vs regrowth")
ax.grid(alpha=.3); ax.legend(loc="upper right", fontsize=9)
ax.text(0.02, 0.02, f"log-rank p={lr.p_value:.1e}\nHR={hr['hr']:.2f} [{hr['lo']:.2f},{hr['hi']:.2f}]",
        transform=ax.transAxes, fontsize=9, va="bottom",
        bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=.8))
fig.tight_layout(); fig.savefig(os.path.join(OUT, "km_moller_taxonomy_regrowth.png")); plt.close()

# ---- save ----
summary = {
    "n_patients": int(len(df)), "n_events": int(df["event"].sum()),
    "event_rate_pct": round(100 * df["event"].mean(), 1),
    "per_k": json.loads(tbl.reset_index().to_json(orient="records")),
    "high_vs_low": {"logrank_p": float(lr.p_value), "cox_hr": hr["hr"],
                    "ci95": [hr["lo"], hr["hi"]], "cox_p": hr["p"]},
    "k4_vs_rest": {"logrank_p": float(lr4.p_value), "cox_hr": hr4["hr"],
                   "ci95": [hr4["lo"], hr4["hi"]], "cox_p": hr4["p"]},
    "regrow_pct_kruskal_p": kw.pvalue,
}
with open(os.path.join(OUT, "prognostic_validation_summary.json"), "w") as f:
    json.dump(summary, f, indent=2)
df.to_csv(os.path.join(OUT, "moller_patient_timeevent.csv"), index=False)
print("\nSaved ->", OUT)
print(json.dumps(summary, indent=2))
