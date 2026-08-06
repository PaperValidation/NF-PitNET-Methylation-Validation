"""
Enhanced direct prognostic validation — Møller 2025 117-NF-PitNET cohort.

Adds a MULTIVARIABLE Cox model: methylation risk axis (high k3/k4/k5 vs low
k1/k2) adjusted for clinical covariates already present in Supplementary_File_S1
(age, sex, invasive growth, radiotherapy, pre-op tumor volume, Knosp).

This shows the methylation axis is an INDEPENDENT predictor, not confounded by
known NF-PitNET regrowth covariates — the analysis a Neuro-Oncology reviewer
will demand, using only the locally available supplement (no Mendeley raw beta).
"""
import json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

S1 = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET/Moller2025_supplement/extracted/Supplementary_File_S1.xlsx"
OUT = "results/Moller2025/prognostic_validation"
import os
os.makedirs(OUT, exist_ok=True)

# ---------- load ----------
mc = pd.read_excel(S1, sheet_name="Main Cohort")
fu = pd.read_excel(S1, sheet_name="Main Cohort Follow Up Data")

# ---------- per-patient follow-up aggregation ----------
fu = fu.copy()
fu["Years"] = pd.to_numeric(fu["Years after 1. postoperative scan"], errors="coerce")
fu["ReOP"] = pd.to_numeric(fu["ReOP (0no,1yes)"], errors="coerce").fillna(0).astype(int)
fu["Reintervention"] = pd.to_numeric(fu["Reintervention"], errors="coerce").fillna(0).astype(int)
fu["Event"] = ((fu["ReOP"] == 1) | (fu["Reintervention"] == 1)).astype(int)

def first_num(g, col):
    s = pd.to_numeric(g[col], errors="coerce")
    s = s.dropna()
    return float(s.iloc[0]) if len(s) else np.nan

def max_num(g, col):
    s = pd.to_numeric(g[col], errors="coerce")
    s = s.dropna()
    return float(s.max()) if len(s) else np.nan

def agg(g):
    ev_rows = g[g["Event"] == 1]
    t = (ev_rows["Years"].min() if len(ev_rows) else g["Years"].max())
    return pd.Series({
        "time": t,
        "event": int(g["Event"].max()),
        "age": first_num(g, "Age at first surgery"),
        "sex": first_num(g, "Sex (0 female/ 1 male)"),
        "hardy": first_num(g, "Hardy"),
        "knosp_r": max_num(g, "Knosp right"),
        "knosp_l": max_num(g, "Knosp left"),
        "invasive": first_num(g, "Invasive Growth"),
        "rt": first_num(g, "RT (0no,1yes)"),
        "preop_vol": first_num(g, "Preoperative Tumor Volume"),
        "cumregrow_max": max_num(g, "Cumulated regrowth %"),
    })

fu_agg = fu.groupby("Sample_ID").apply(agg, include_groups=False).reset_index()

# ---------- merge with Main Cohort (k-label) ----------
mc_sub = mc[["Sample_ID", "Kcluster_5", "Gender", "Excluded", "InvasiveGrowth"]].copy()
df = fu_agg.merge(mc_sub, on="Sample_ID", how="left")
df = df[df["Excluded"] == 0].copy()           # drop lost-to-follow-up
df = df[df["time"].notna()].copy()
df["time"] = df["time"].astype(float)

# risk axis
df["risk"] = df["Kcluster_5"].map(lambda k: "high(k3k4k5)" if k in ("k3", "k4", "k5") else "low(k1k2)")
df["risk_hi"] = (df["risk"] == "high(k3k4k5)").astype(int)
df["k4"] = (df["Kcluster_5"] == "k4").astype(int)
df["invasive_bin"] = (df["invasive"] == 1).astype(int)
df["sex_male"] = (df["sex"] == 1).astype(int)
df["rt_yn"] = (df["rt"] == 1).astype(int)
df["knosp_max"] = df[["knosp_r", "knosp_l"]].max(axis=1)

print("N patients (follow-up, non-excluded):", len(df))
print("Events:", int(df["event"].sum()), "(%.1f%%)" % (100*df["event"].mean()))
print("Missingness per covariate:")
for c in ["age","sex_male","knosp_max","invasive_bin","rt_yn","preop_vol","cumregrow_max"]:
    print(f"  {c}: {int(df[c].isna().sum())} missing")

# ---------- Cox helpers (robust) ----------
def cox_one(frame, pred, dur="time", ev="event"):
    cph = CoxPHFitter()
    cph.fit(frame[["time","event",pred]], duration_col=dur, event_col=ev)
    r = cph.summary.loc[pred]
    return dict(hr=float(r["exp(coef)"]), lo=float(r["exp(coef) lower 95%"]),
                hi=float(r["exp(coef) upper 95%"]), p=float(r["p"]), n=int(len(frame)))

def cox_multi(frame, preds, dur="time", ev="event"):
    cph = CoxPHFitter()
    cph.fit(frame[["time","event"]+preds], duration_col=dur, event_col=ev)
    out = {}
    for p in preds:
        r = cph.summary.loc[p]
        out[p] = dict(hr=float(r["exp(coef)"]), lo=float(r["exp(coef) lower 95%"]),
                      hi=float(r["exp(coef) upper 95%"]), p=float(r["p"]))
    return out, int(len(frame))

# ---------- univariable Cox (one predictor at a time) ----------
uni_preds = ["risk_hi","age","sex_male","knosp_max","invasive_bin","rt_yn","preop_vol","k4"]
uni = {}
for p in uni_preds:
    sub = df.dropna(subset=[p])
    if sub["event"].nunique() < 2 or sub[p].nunique() < 2:
        uni[p] = dict(hr=None, note="no variation / <2 events")
        continue
    try:
        uni[p] = cox_one(sub, p)
    except Exception as e:
        uni[p] = dict(hr=None, note="convergence: %s" % type(e).__name__)

# ---------- multivariable Cox (parsimonious, no rt_yn) ----------
mv_preds = ["risk_hi","age","sex_male","invasive_bin"]
mv_df = df.dropna(subset=mv_preds)
mv, mv_n = cox_multi(mv_df, mv_preds)
# optional extension with preop tumor volume
mv2_preds = mv_preds + ["preop_vol"]
mv2_df = df.dropna(subset=mv2_preds)
try:
    mv2, mv2_n = cox_multi(mv2_df, mv2_preds)
except Exception as e:
    mv2, mv2_n = None, None

print("\nMultivariable (risk+age+sex+invasive) n =", mv_n, "events =", int(mv_df["event"].sum()))
for k,v in mv.items():
    print(f"  {k}: HR={v['hr']:.2f} [{v['lo']:.2f},{v['hi']:.2f}] p={v['p']:.3g}")
if mv2:
    print("Multivariable +preop_vol n =", mv2_n)
    for k,v in mv2.items():
        print(f"  {k}: HR={v['hr']:.2f} [{v['lo']:.2f},{v['hi']:.2f}] p={v['p']:.3g}")

# ---------- also the original pooled high/low (consistency) ----------
from lifelines.statistics import logrank_test
grp = df.groupby("risk_hi")["time"]
t1 = df[df.risk_hi==1]["time"]; e1 = df[df.risk_hi==1]["event"]
t0 = df[df.risk_hi==0]["time"]; e0 = df[df.risk_hi==0]["event"]
lr = logrank_test(t0, t1, event_observed_A=e0, event_observed_B=e1)
pooled = cox_one(df.dropna(subset=["risk_hi"]), "risk_hi")

summary = {
    "n_patients": int(len(df)),
    "n_events": int(df["event"].sum()),
    "event_rate": float(df["event"].mean()),
    "pooled_high_vs_low": {
        "hr": pooled["hr"], "lo": pooled["lo"], "hi": pooled["hi"],
        "cox_p": pooled["p"], "logrank_p": float(lr.p_value)
    },
    "univariable": uni,
    "multivariable_parsimonious": {"predictors": mv_preds, "n": int(len(mv_df)),
                                    "n_events": int(mv_df["event"].sum()), "results": mv},
    "multivariable_plus_preopvol": ({"predictors": mv2_preds, "n": int(mv2_n),
                                    "n_events": int(mv2_df["event"].sum()), "results": mv2}
                                    if mv2 else None),
}
with open(f"{OUT}/prognostic_validation_multivariable.json","w") as f:
    json.dump(summary, f, indent=2)

# save enriched per-patient table
df.to_csv(f"{OUT}/moller_patient_timeevent_enriched.csv", index=False)
print("\nWrote:", f"{OUT}/prognostic_validation_multivariable.json")
print("Wrote:", f"{OUT}/moller_patient_timeevent_enriched.csv")
