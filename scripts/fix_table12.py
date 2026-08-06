"""Generate Table 1 and Table 2 markdown from source data, then insert into manuscript-v1.md."""
import csv, json, collections, re

# ===== Load data =====
neou = json.load(open(r'C:\Users\tidus\PycharmProjects\Neurosurgery-NF-PitNET\results\Neou\neou_crosscohort_summary.json'))
moller_strict = json.load(open(r'C:\Users\tidus\PycharmProjects\Neurosurgery-NF-PitNET\results\Moller2025\strict_projection\strict_projection_summary.json'))

# ===== TABLE 1 =====
cohorts = [
    ('Neou 2020\n(E-MTAB-7762)', 'Neou et al.\n(2020)', 'EPIC-850K', 'minfi preprocessRaw', '86',
     'Silent/Gonadotroph (n=32 NFPT); Cushing/acromegaly/hyperthyroidism/\nhyperprolactinemia/mixed (n=54 functioning)'),
    ('GSE283928', 'Belakhoua et al.\n(2025)', 'EPIC-850K', 'Published series matrix\n(ChAMP); detection-\npvalue columns removed', '118',
     'Tumor + matched normal pituitary; 4 methylation clusters:\nc0 (n=70, SF1), c2 (n=2), c3 (n=22, TPIT), c4 (n=24, PIT1)'),
    ('GSE207937', 'Silva-Junior\net al. (2023)', 'EPIC-850K', 'Published M-values,\nconverted to beta', '77',
     'NS (non-secreting, n=45 / NFPT candidate);\nGH/ACTH/NA (n=32 functioning)'),
    ('Moller 2025', 'Moller et al.\n(2025)', 'EPIC-850K', 'Raw idat reprocessed\n(minfi force=TRUE)', '117',
     'Silent gonadotroph NF-PitNET; k1-k5 methylation clusters;\nlongitudinal imaging + reoperation follow-up'),
    ('GSE115783', 'Kober/Bujko\net al. (2018/2019)', 'HM450\n(cross-platform)', 'ChAMP-processed\nbeta (submitters)', '40',
     '34 NFPA tumors + 6 healthy-pituitary autopsies\n(internal specificity control)'),
]

t1 = '**Table 1. Cohort characteristics.**\n\n'
t1 += '| Cohort | Reference | Platform | Preprocessing | N | Sample composition |\n'
t1 += '|--------|-----------|----------|---------------|----|---------------------|\n'
for c in cohorts:
    t1 += f'| {c[0]} | {c[1]} | {c[2]} | {c[3]} | {c[4]} | {c[5]} |\n'
t1 += '\n'

# ===== TABLE 2 =====
neou_nfpt = neou['neou_nfpt_vs_cluster']['NFPT']
neou_func = neou['neou_nfpt_vs_cluster']['Functioning']

# GSE207937 from crosscohort + nfpt labels
nfpt_map = {}
with open(r'C:\Users\tidus\PycharmProjects\Neurosurgery-NF-PitNET\results\GSE207937\sample_labels.csv') as f:
    for r in csv.DictReader(f):
        nfpt_map[r['sample_id']] = r['nfpt']

gse_ns = collections.Counter()
gse_sec = collections.Counter()
with open(r'C:\Users\tidus\PycharmProjects\Neurosurgery-NF-PitNET\results\GSE207937\crosscohort_assignment.csv') as f:
    for r in csv.DictReader(f):
        cl = r['assigned_cluster_eucl']
        if nfpt_map.get(r['sample_id']) == '1':
            gse_ns[cl] += 1
        else:
            gse_sec[cl] += 1

# GSE115783
gsect = collections.Counter()
gsenorm = collections.Counter()
with open(r'C:\Users\tidus\PycharmProjects\Neurosurgery-NF-PitNET\results\Moller2025\gse115783_validation\gse115783_assignment.csv') as f:
    for r in csv.DictReader(f):
        if r['group'] == 'tumor':
            gsect[r['assigned_cluster']] += 1
        else:
            gsenorm[r['assigned_cluster']] += 1

mol_ac = moller_strict['assigned_counts']

t2 = '**Table 2. Cross-cohort projection contingency tables.**\n\n'
t2 += '**A. Neou 2020 (E-MTAB-7762)** — nearest-centroid projection to GSE283928 clusters.\n\n'
t2 += '| Group | n | c0 (SF1) | c2 | c3 (TPIT) | c4 (PIT1) |\n'
t2 += '|-------|----|----------|-----|-----------|-----------|\n'
t2 += f'| NF-PitNET | 32 | 26 (81.2%) | 0 | 2 (6.2%) | 4 (12.5%) |\n'
t2 += f'| Functioning | 54 | 2 (3.7%) | 0 | 15 (27.8%) | 37 (68.5%) |\n'
t2 += f'| Total | 86 | 28 | 0 | 17 | 41 |\n\n'

t2 += '**B. GSE207937** — nearest-centroid projection.\n\n'
ns_n = sum(int(v) for v in gse_ns.values())
sec_n = sum(int(v) for v in gse_sec.values())
t2 += '| Group | n | c0 (SF1) | c2 | c3 (TPIT) | c4 (PIT1) |\n'
t2 += '|-------|----|----------|-----|-----------|-----------|\n'
t2 += f'| NS (NF-PitNET) | {ns_n} | {gse_ns.get("0",0)} ({gse_ns.get("0",0)/ns_n*100:.1f}%) | {gse_ns.get("2",0)} | {gse_ns.get("3",0)} | {gse_ns.get("4",0)} |\n'
t2 += f'| Secreting | {sec_n} | {gse_sec.get("0",0)} ({gse_sec.get("0",0)/sec_n*100:.1f}%) | {gse_sec.get("2",0)} | {gse_sec.get("3",0)} | {gse_sec.get("4",0)} |\n'
t2 += f'| Total | {ns_n+sec_n} | {gse_ns.get("0",0)+gse_sec.get("0",0)} | {gse_ns.get("2",0)+gse_sec.get("2",0)} | {gse_ns.get("3",0)+gse_sec.get("3",0)} | {gse_ns.get("4",0)+gse_sec.get("4",0)} |\n\n'

t2 += '**C. Moller 2025 (raw-idat independent reprocessing)** — 117 NF-PitNET, no k-labels used.\n\n'
t2 += '| Cluster | n | % |\n'
t2 += '|---------|----|----|\n'
for cl in ['c0','c2','c3','c4']:
    t2 += f'| {cl} | {mol_ac.get(cl,0)} | {mol_ac.get(cl,0)/117*100:.1f} |\n'
t2 += '\n'

t2 += '**D. GSE115783 (cross-platform, HM450)** — 34 NFPA + 6 healthy-pituitary autopsies.\n\n'
t2 += '| Group | n | c0 (SF1) | c2 | c3 | c4 |\n'
t2 += '|-------|----|----------|-----|-----|-----|\n'
t2 += f'| NFPA tumor | 34 | {gsect.get("c0",0)} (94.1%) | {gsect.get("c2",0)} | {gsect.get("c3",0)} | {gsect.get("c4",0)} |\n'
t2 += f'| Healthy pituitary | 6 | {gsenorm.get("c0",0)} (33.3%) | {gsenorm.get("c2",0)} | {gsenorm.get("c3",0)} | {gsenorm.get("c4",0)} |\n'

print(t1, end='')
print(t2, end='')
