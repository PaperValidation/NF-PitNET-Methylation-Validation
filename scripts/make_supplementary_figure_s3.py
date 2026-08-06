#!/usr/bin/env python
# Supplementary Figure S3: Drug-repositioning candidates from Enrichr/LINCS
# Regenerates the CMap connectivity-mapping figure with explicit
# "Supplementary Figure S3" labeling and evidence-tier coloring on the
# drug-class frequency panel (honest grading carried into the figure).
import re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

OUT = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/results/CMap"
rank = pd.read_csv(f"{OUT}/drug_reversal_ranking.csv")
top = pd.read_csv(f"{OUT}/drug_reversal_top40.csv")


def clean(c):
    c = str(c).replace('6H-', '').replace('24H-', '').replace('3H-', '').replace('12H-', '').replace('48H-', '')
    c = re.sub(r'-\d+(\.\d+)?$', '', c)
    return c.strip('-')


CLASS_RULES = {
    'BET/BRD': ['BRD-', 'JQ1', 'OTX', 'CPI-', 'GSK-', 'ABBV-'],
    'MEK': ['selumetinib', 'PD-0325901', 'trametinib', 'binimetinib', 'cobimetinib', 'PD0325901', 'AZD-8330', 'pimasertib'],
    'mTOR/PI3K': ['AZD-8055', 'everolimus', 'temsirolimus', 'rapamycin', 'torin', 'PF-', 'GSK-2126458', 'BEZ235', 'buparlisib', 'alpelisib'],
    'SFK/ABL': ['dasatinib', 'nilotinib', 'bosutinib', 'imatinib', 'saracatinib'],
    'HDAC': ['vorinostat', 'valproic', 'trichostatin', 'romidepsin', 'panobinostat', 'entinostat', 'belinostat'],
    'EGFR/HER': ['erlotinib', 'gefitinib', 'lapatinib', 'SU-11274', 'foretinib', 'afatinib', 'vandetanib', 'cetuximab'],
    'RAF/MEK': ['vemurafenib', 'dabrafenib', 'zelboraf', 'encorafenib'],
    'CDK': ['palbociclib', 'dinaciclib', 'roscovitine', 'seliciclib', 'CDK'],
    'JAK/STAT': ['ruxolitinib', 'tofacitinib', 'fedratinib', 'momelotinib'],
    'HSP90': ['radicicol', 'geldanamycin', '17-AAG', 'tanespimycin', 'AUY'],
    'Glucocorticoid': ['dexamethasone', 'prednisolone', 'hydrocortisone'],
    'Hormone': ['estradiol', 'progesterone', 'letrozole', 'tamoxifen', 'diethylstilbestrol'],
}


def classify(comp):
    c = clean(comp)
    for cls, kws in CLASS_RULES.items():
        for kw in kws:
            if kw.lower() in c.lower():
                return cls
    return 'Other'


# Evidence-tier color map (PitNET translational evidence, from deep-research audit)
TIER_COLOR = {
    'HDAC': '#1b7837',          # deepest
    'mTOR/PI3K': '#1b7837',     # deepest
    'BET/BRD': '#2a9d8f',       # strong NFPA models
    'MEK': '#e76f51',           # corticotroph-only (weak)
    'RAF/MEK': '#e76f51',       # MEK-pathway, corticotroph-limited
    'CDK': '#c1121f',           # single NFPA case (lowest)
    'SFK/ABL': '#8d99ae',
    'EGFR/HER': '#8d99ae',
    'HSP90': '#8d99ae',
    'JAK/STAT': '#8d99ae',
    'Glucocorticoid': '#8d99ae',
    'Hormone': '#8d99ae',
    'Other': '#8d99ae',
}
TIER_LABEL = {
    'HDAC': 'Deepest PitNET evidence',
    'mTOR/PI3K': 'Deepest PitNET evidence',
    'BET/BRD': 'Strong NFPA models',
    'MEK': 'Corticotroph-only (weak)',
    'RAF/MEK': 'Corticotroph-only (weak)',
    'CDK': 'Single NFPA case (lowest)',
    'Other': 'Other / non-prioritized',
}

rank['compound_clean'] = rank['compound'].map(clean)
rank['class'] = rank['compound_clean'].map(classify)
cls_counts = rank['class'].value_counts()
cls_counts = cls_counts[cls_counts.index != 'Other'].head(12)

# Panel A: top 20 individual instances
ta = top.copy()
ta['label'] = ta['compound'].map(clean)
ta = ta.sort_values('reversal_score').tail(20)
colorsA = [('#2a9d8f' if r.rev_up_cs >= r.rev_down_cs else '#e76f51') for _, r in ta.iterrows()]

fig, (axA, axB) = plt.subplots(1, 2, figsize=(15, 7.8))
fig.suptitle('Supplementary Figure S3. Drug-repositioning candidates from Enrichr/LINCS connectivity mapping\n(exploratory; NFPT transcriptional signature vs LINCS L1000 perturbation gene sets)',
             fontsize=12, fontweight='bold')

axA.barh(ta['label'], ta['reversal_score'], color=colorsA)
axA.set_title('A. Top 20 reversal candidates (individual LINCS L1000 instances)', fontsize=11, fontweight='bold')
axA.set_xlabel('Reversal score (Σ reversal − Σ pro-disease; Enrichr Combined Score)')
axA.tick_params(axis='y', labelsize=8)
axA.legend(handles=[Patch(color='#2a9d8f', label='driven by down-regulating NFPT-UP genes'),
                    Patch(color='#e76f51', label='driven by up-regulating NFPT-DOWN genes')],
            fontsize=8, loc='lower right')

# Panel B: drug-class frequency, colored by evidence tier
bar_colors = [TIER_COLOR.get(c, '#8d99ae') for c in cls_counts.index[::-1]]
axB.barh(cls_counts.index[::-1], cls_counts.values[::-1], color=bar_colors)
axB.set_title('B. Drug-class frequency among reversal-positive perturbations\n(n=%d instances)' % len(rank), fontsize=11, fontweight='bold')
axB.set_xlabel('Number of LINCS instances reversing the NFPT program')
axB.tick_params(axis='y', labelsize=9)
tier_handles = [
    Patch(color='#1b7837', label='HDAC / mTOR — deepest PitNET evidence'),
    Patch(color='#2a9d8f', label='BET/BRD — strong NFPA models'),
    Patch(color='#e76f51', label='MEK / RAF-MEK — corticotroph-only (weak)'),
    Patch(color='#c1121f', label='CDK4/6 — single NFPA case (lowest)'),
    Patch(color='#8d99ae', label='Other / non-prioritized classes'),
]
axB.legend(handles=tier_handles, fontsize=7.5, loc='lower right')

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(f"{OUT}/supplementary_figure_s3.png", dpi=140, bbox_inches='tight')
print("saved", f"{OUT}/supplementary_figure_s3.png")
print("class_counts:")
print(cls_counts.to_string())
