"""Forest plot of Cox HRs from the Møller 117 multivariable validation."""
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

J = "results/Moller2025/prognostic_validation/prognostic_validation_multivariable.json"
OUT = "results/Moller2025/prognostic_validation/forest_moller_cox.png"
s = json.load(open(J))

rows = []  # (label, hr, lo, hi, significant)
def add(label, d, sig_thr=0.05):
    if d and d.get("hr") is not None:
        rows.append((label, d["hr"], d["lo"], d["hi"], d["p"] < sig_thr))

add("Methylation risk axis (uni)", s["univariable"]["risk_hi"])
add("Methylation risk axis (adj)", s["multivariable_parsimonious"]["results"]["risk_hi"])
add("Methylation risk axis (+vol)", s["multivariable_plus_preopvol"]["results"]["risk_hi"])
add("Invasive growth (uni)", s["univariable"]["invasive_bin"])
add("Invasive growth (adj)", s["multivariable_parsimonious"]["results"]["invasive_bin"])
add("Pre-op tumor volume (uni)", s["univariable"]["preop_vol"])
add("Pre-op tumor volume (+vol)", s["multivariable_plus_preopvol"]["results"]["preop_vol"])
add("Age at surgery (uni)", s["univariable"]["age"])
add("Male sex (uni)", s["univariable"]["sex_male"])
add("Knosp max (uni)", s["univariable"]["knosp_max"])
add("k4 vs rest (uni)", s["univariable"]["k4"])

labels = [r[0] for r in rows]
hrs = [r[1] for r in rows]
los = [r[2] for r in rows]
his = [r[3] for r in rows]
sig = [r[4] for r in rows]

y = np.arange(len(rows))[::-1]
fig, ax = plt.subplots(figsize=(7.2, 5.0))
for i, (hr, lo, hi, s_) in enumerate(zip(hrs, los, his, sig)):
    color = "#c0392b" if s_ else "#555555"
    ax.plot([lo, hi], [y[i], y[i]], color=color, lw=1.8, zorder=1)
    ax.plot(hr, y[i], "o", color=color, ms=6, zorder=2)
ax.axvline(1.0, color="black", lw=1, ls="--")
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel("Hazard ratio (95% CI)", fontsize=10)
ax.set_xscale("log")
ax.set_xlim(0.4, 6)
ax.set_title("Cox regression — regrowth (reoperation/reintervention)\nMøller 117 NF-PitNET cohort", fontsize=10)
# annotate HR values
for i, (hr, lo, hi, s_) in enumerate(zip(hrs, los, his, sig)):
    tag = f"{hr:.2f} [{lo:.2f}, {hi:.2f}]" + ("*" if s_ else "")
    ax.text(hi*1.05, y[i], tag, va="center", fontsize=7.5, color="#333333")
ax.set_xlim(0.4, 9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(OUT, dpi=200)
print("Wrote", OUT)
