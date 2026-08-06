#!/usr/bin/env python
# Figure: CMap/Enrichr drug-repositioning results for NFPT signature
import json, re
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT="C:/Users/tidus/PycharmProjects/Neurosurgery-NF-PitNET/results/CMap"
rank=pd.read_csv(f"{OUT}/drug_reversal_ranking.csv")
top=pd.read_csv(f"{OUT}/drug_reversal_top40.csv")

def clean(c):
    c=str(c).replace('6H-','').replace('24H-','').replace('3H-','').replace('12H-','').replace('48H-','')
    c=re.sub(r'-\d+(\.\d+)?$','',c)
    return c.strip('-')

CLASS_RULES={
 'BET/BRD':['BRD-','JQ1','OTX','CPI-','GSK-','ABBV-'],
 'MEK':['selumetinib','PD-0325901','trametinib','binimetinib','cobimetinib','PD0325901','AZD-8330','pimasertib'],
 'mTOR/PI3K':['AZD-8055','everolimus','temsirolimus','rapamycin','torin','PF-','GSK-2126458','BEZ235','buparlisib','alpelisib'],
 'SFK/ABL':['dasatinib','nilotinib','bosutinib','imatinib','saracatinib','dasatinib'],
 'HDAC':['vorinostat','valproic','trichostatin','romidepsin','panobinostat','entinostat','belinostat'],
 'EGFR/HER':['erlotinib','gefitinib','lapatinib','SU-11274','foretinib','afatinib','vandetanib','cetuximab'],
 'RAF/MEK':['vemurafenib','dabrafenib','zelboraf','encorafenib'],
 'CDK':['palbociclib','dinaciclib','roscovitine','seliciclib','CDK'],
 'JAK/STAT':['ruxolitinib','tofacitinib','fedratinib','momelotinib'],
 'HSP90':['radicicol','geldanamycin','17-AAG','tanespimycin','AUY'],
 'Glucocorticoid':['dexamethasone','prednisolone','hydrocortisone'],
 'Hormone':['estradiol','progesterone','letrozole','tamoxifen','diethylstilbestrol'],
}
def classify(comp):
    c=clean(comp)
    for cls,kws in CLASS_RULES.items():
        for kw in kws:
            if kw.lower() in c.lower(): return cls
    return 'Other'
rank['compound_clean']=rank['compound'].map(clean)
rank['class']=rank['compound_clean'].map(classify)
# class frequency over all positive-reversal instances
cls_counts=rank['class'].value_counts()
cls_counts=cls_counts[cls_counts.index!='Other'].head(12)

# Panel A: top 20 individual instances
ta=top.copy()
ta['label']=ta['compound'].map(clean)
ta=ta.sort_values('reversal_score').tail(20)
# color by dominant direction
colors=[('#2a9d8f' if r.rev_up_cs>=r.rev_down_cs else '#e76f51') for _,r in ta.iterrows()]

fig,(axA,axB)=plt.subplots(1,2,figsize=(15,7.5))
axA.barh(ta['label'],ta['reversal_score'],color=colors)
axA.set_title('Top 20 reversal candidates (NFPT signature → LINCS L1000)',fontsize=12,fontweight='bold')
axA.set_xlabel('Reversal score (Σ reversal − Σ pro-disease, Enrichr Combined Score)')
axA.tick_params(axis='y',labelsize=8)
from matplotlib.patches import Patch
axA.legend(handles=[Patch(color='#2a9d8f',label='driven by down-regulating NFPT-UP genes'),
                    Patch(color='#e76f51',label='driven by up-regulating NFPT-DOWN genes')],
           fontsize=8,loc='lower right')

axB.barh(cls_counts.index[::-1],cls_counts.values[::-1],color='#264653')
axB.set_title('Drug-class frequency among all reversal-positive\nperturbations (n=%d instances)'%len(rank),fontsize=12,fontweight='bold')
axB.set_xlabel('Number of LINCS instances reversing NFPT program')
axB.tick_params(axis='y',labelsize=9)

plt.tight_layout()
plt.savefig(f"{OUT}/cmap_figure.png",dpi=130,bbox_inches='tight')
print("saved",f"{OUT}/cmap_figure.png")
print("\nTop drug classes (reversal-positive instances):")
print(cls_counts.to_string())
# also report MEK/mTOR/BET convergence note
for cls in ['MEK','mTOR/PI3K','BET/BRD','SFK/ABL','HDAC']:
    n=int((rank['class']==cls).sum())
    print(f"  {cls}: {n} instances")
json.dump({'class_counts':cls_counts.to_dict(),
           'n_positive':int(len(rank)),
           'top20':ta[['instance','label','reversal_score']].to_dict('records')},
          open(f"{OUT}/cmap_summary.json","w"),indent=1)
