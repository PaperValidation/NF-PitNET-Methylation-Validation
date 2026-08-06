#!/usr/bin/env python
# CMap-style connectivity / drug-repositioning via Enrichr LINCS-L1000 + CMap drug-perturbation libraries.
# NFPT disease signature (GSE209903 FPT_vs_NFPT, DESeq2) queried against drug perturbation signatures.
# Reversal logic: a drug that DOWN-regulates genes up in NFPT (and UP-regulates genes down in NFPT)
# is a candidate that reverses the NFPT transcriptional program.
import os, json, re
os.environ['HTTPS_PROXY'] = 'socks5h://127.0.0.1:8428'
os.environ['HTTP_PROXY']  = 'socks5h://127.0.0.1:8428'

import gseapy as gp
import pandas as pd

OUT = "C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/results/CMap"
sym = json.load(open(f"{OUT}/nfpt_symbols.json"))
up   = sym['up_strict']      # genes up in NFPT (disease-up)
down = sym['down_strict']    # genes down in NFPT (disease-down)
print(f"signature: up={len(up)} down={len(down)}")

def compound_of(term):
    parts = term.split()
    cands = [p for p in parts if '-' in p and any(c.isdigit() for c in p) and any(c.isalpha() for c in p)]
    return cands[-1] if cands else parts[-1]

def run(query_list, library, tag):
    try:
        e = gp.enrichr(gene_list=query_list, gene_sets=library,
                       organism='human', outdir=None, cutoff=0.05)
        df = e.results
        df['_tag'] = tag
        df['_compound'] = df['Term'].map(compound_of)
        df.to_csv(f"{OUT}/enrichr_{tag}.csv", index=False)
        print(f"  [{tag}] {library}: {len(df)} significant sets")
        return df
    except Exception as ex:
        print(f"  [{tag}] {library} FAILED: {str(ex)[:160]}")
        return pd.DataFrame()

# Directional LINCS L1000 chemical perturbations
res = {}
res['up_vs_chem_down'] = run(up,   'LINCS_L1000_Chem_Pert_down', 'up_vs_chem_down')  # drug downregulates disease-UP -> reversal (+)
res['up_vs_chem_up']   = run(up,   'LINCS_L1000_Chem_Pert_up',   'up_vs_chem_up')    # drug upregulates disease-UP -> pro-disease (-)
res['down_vs_chem_up'] = run(down, 'LINCS_L1000_Chem_Pert_up',   'down_vs_chem_up')  # drug upregulates disease-DOWN -> reversal (+)
res['down_vs_chem_down']=run(down, 'LINCS_L1000_Chem_Pert_down', 'down_vs_chem_down')# drug downregulates disease-DOWN -> pro-disease (-)
# Classic CMap-style (GEO drug perturbations, directional-agnostic overlap)
res['up_vs_cmap']      = run(up,   'Drug_Perturbations_from_GEO_2014', 'up_vs_cmap')
res['down_vs_cmap']    = run(down, 'Drug_Perturbations_from_GEO_2014', 'down_vs_cmap')

def cs_map(df):
    if df.empty: return {}
    out = {}
    for _,r in df[['Term','_compound','Combined Score','Adjusted P-value','Odds Ratio','Genes']].iterrows():
        out[r['Term']] = r.to_dict()
    return out

m = {k: cs_map(v) for k,v in res.items()}

drugs = set()
for k in ['up_vs_chem_down','up_vs_chem_up','down_vs_chem_up','down_vs_chem_down']:
    drugs |= set(m[k].keys())
print(f"\nTotal directional LINCS instances evaluated: {len(drugs)}")

def get(d, k, key, default=0):
    v = m[k].get(d, {})
    x = v.get(key, default)
    return x if pd.notna(x) else default

rows=[]
for d in drugs:
    rev_up    = get(d,'up_vs_chem_down','Combined Score',0) or 0   # drug downregulates disease-UP genes
    prod_up   = get(d,'up_vs_chem_up',  'Combined Score',0) or 0   # drug upregulates disease-UP genes (bad)
    rev_down  = get(d,'down_vs_chem_up','Combined Score',0) or 0   # drug upregulates disease-DOWN genes
    prod_down = get(d,'down_vs_chem_down','Combined Score',0) or 0 # drug downregulates disease-DOWN genes (bad)
    reversal = rev_up + rev_down - prod_up - prod_down
    if reversal <= 0: continue
    rows.append({'instance':d,'compound':m['up_vs_chem_down'].get(d,{}).get('_compound') or m['down_vs_chem_up'].get(d,{}).get('_compound') or m['up_vs_chem_up'].get(d,{}).get('_compound') or m['down_vs_chem_down'].get(d,{}).get('_compound'),
                 'reversal_score':round(reversal,2),
                 'rev_up_cs':round(rev_up,1),'prod_up_cs':round(prod_up,1),
                 'rev_down_cs':round(rev_down,1),'prod_down_cs':round(prod_down,1)})
rank = pd.DataFrame(rows).sort_values('reversal_score', ascending=False)
rank.to_csv(f"{OUT}/drug_reversal_ranking.csv", index=False)
rank.head(40).to_csv(f"{OUT}/drug_reversal_top40.csv", index=False)
print(f"Drugs with positive reversal score: {len(rank)}")
print("\n=== TOP 25 REVERSAL CANDIDATES (LINCS directional) ===")
print(rank.head(25).to_string(index=False))

# CMap-style overlap (directional-agnostic): disease-up/down genes overlap drug perturbation sets
cmap_rows=[]
for tag,df in [('up_vs_cmap',res['up_vs_cmap']),('down_vs_cmap',res['down_vs_cmap'])]:
    if df.empty: continue
    for _,r in df.iterrows():
        genes = str(r['Genes']).split(';') if pd.notna(r['Genes']) else []
        cmap_rows.append({'drug':r['Term'],'direction':tag,'adj_p':r['Adjusted P-value'],
                          'odds':r['Odds Ratio'],'combined_score':r['Combined Score'],'n_overlap':len(genes)})
cmap_df=pd.DataFrame(cmap_rows).sort_values('combined_score',ascending=False)
cmap_df.to_csv(f"{OUT}/cmap_style_overlap.csv", index=False)
print(f"\n=== TOP 20 CMap-style (disease-up/down overlap) ===")
print(cmap_df.head(20).to_string(index=False))
print("\nDONE")
