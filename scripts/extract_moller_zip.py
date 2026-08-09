#!/usr/bin/env python3
# One-off: extract Møller 2025 Mendeley zip (flatten top folder) so that
#   Mendeley_Moller2025/Deidentified_IDATs/<Sentrix_ID>/...idat
#   Mendeley_Moller2025/SampleSheet.xlsx
# exist directly (what moller_idat_to_beta.R expects).
import os, zipfile, sys

ZIP = r"G:\CMCC Sync\Medical Data\Neurosurgery\NF-PitNET\Mendeley_Moller2025\Raw DNA Methylation Data; Møller et al., 2025.zip"
DST = r"G:\CMCC Sync\Medical Data\Neurosurgery\NF-PitNET\Mendeley_Moller2025"

with zipfile.ZipFile(ZIP) as zf:
    names = zf.namelist()
    top = names[0].split('/')[0]
    prefix = top + '/'
    targets = [n for n in names if n.startswith(prefix) and not n.endswith('/')]
    total = len(targets)
    print(f"top folder: {top!r}  files to extract: {total}")
    for i, n in enumerate(targets, 1):
        rel = n[len(prefix):]                      # strip top folder
        if not rel:
            continue
        out = os.path.join(DST, rel)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with zf.open(n) as src, open(out, 'wb') as fh:
            fh.write(src.read())
        if i % 50 == 0 or i == total:
            print(f"  {i}/{total}  -> {rel}")
print("EXTRACTION DONE")
# verify
idat = sum(len(fs) for _, _, fs in os.walk(os.path.join(DST, 'Deidentified_IDATs')))
print("Deidentified_IDATs idat count:", idat)
print("SampleSheet exists:", os.path.exists(os.path.join(DST, 'SampleSheet.xlsx')))
