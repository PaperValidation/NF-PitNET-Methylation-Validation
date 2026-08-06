import re, os, json
html = r"G:\CMCC Sync\Medical Data\Neurosurgery\NF-PitNET\pmc13070495_printable.html"
data = open(html, encoding="utf-8").read()
text = re.sub(r"<[^>]+>", " ", data)
text = re.sub(r"&[a-z]+;", " ", text)
text = re.sub(r"\s+", " ", text)

# 1) Extract the cluster-specific promoter gene mentions context
for gene in ["CCDC141","ATP1A3","ANK2","CASC15","HDC","NCKAP5","GUCA1C","RP11-497G19.3","ATP2C1","RSPH6A"]:
    m = re.search(r".{120}"+re.escape(gene)+r".{120}", text)
    print(f"[{gene}]", m.group(0) if m else "not found in body text")

print("\n=== Neou / GSE references ===")
# find sentences mentioning Neou
for m in re.finditer(r".{80}Neou.{120}", text):
    print("-", m.group(0))
# any GSE inside
print("GSEs:", sorted(set(re.findall(r"GSE\d{5,7}", data))))

# 2) Look for supplementary file S2 link / table of DMPs
print("\n=== supp file S2 context ===")
m = re.search(r".{200}Supplementary File S2.{200}", text)
print(m.group(0) if m else "no S2 context")
# try to find any table-like rows with cg probes
cgs = re.findall(r"cg\d{6,8}", data)
print("cg probe ids found in html:", len(cgs), "sample:", cgs[:10])

# 3) Save a structured summary
summary = {
  "paper": "Moller et al. Neuro-Oncology 2025 (noaf269) PMID 41369112",
  "cohorts_external": {"Silva-Junior": "GSE207937 (77, 46 NFPT) + GSE209903 RNAseq",
                        "Belakhoua": "GSE283928 (118, public processed matrix, 12 recurrences)",
                        "Neou": "Cancer Cell 2020 (GSE TBD)"},
  "key_genes": ["CCDC141","ATP1A3","ANK2","CASC15","HDC","NCKAP5","GUCA1C","RP11-497G19.3","ATP2C1","RSPH6A"],
  "n_dmp": 562,
  "clusters": ["k1","k2","k3","k4","k5"],
  "gse_in_article": sorted(set(re.findall(r"GSE\d{5,7}", data)))
}
out = r"G:\CMCC Sync\Medical Data\Neurosurgery\NF-PitNET\noaf269_summary.json"
open(out,"w",encoding="utf-8").write(json.dumps(summary, ensure_ascii=False, indent=2))
print("\nWROTE", out)
