"""Parse GEO series_matrix.txt.gz metadata block into a per-sample table.
Only the '!' header section is read (metadata), not the giant expression table.
Outputs <DATASET>/samples/<GSE>_samples.csv
"""
import os, re, csv, gzip, sys

GDRIVE = "G:/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
OUTDIR = os.path.join(GDRIVE, "samples")
os.makedirs(OUTDIR, exist_ok=True)

# metadata fields we care about (collected in sample order)
WANT = [
    "geo_accession", "title", "source_name_ch1", "organism_ch1",
    "molecule_ch1", "characteristics_ch1", "treatment_protocol_ch1",
    "growth_protocol_ch1", "extract_protocol_ch1", "label_ch1",
    "label_protocol_ch1", "hyb_protocol_ch1", "scan_protocol_ch1",
    "description", "data_processing", "platform_id", "contact_*, "
    "supplementary_file", "series_id", "status", "submission_date",
    "last_update_date", "type",
]

def field_key(line):
    # "!Sample_characteristics_ch1" -> "characteristics_ch1"
    m = re.match(r"^!Sample_([A-Za-z0-9_]+)", line)
    if m:
        return "sample_" + m.group(1)
    m = re.match(r"^!Series_([A-Za-z0-9_]+)", line)
    if m:
        return "series_" + m.group(1)
    return None

def parse(path):
    meta = {}  # key -> list of values per sample
    order = []
    n_samples = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                continue
            key = field_key(line.rstrip("\n"))
            if not key:
                continue
            # first token after '!' is the field name; rest is tab-separated values
            parts = line.rstrip("\n").split("\t")
            # parts[0] = "!Sample_xxx" ; parts[1:] = per-sample values
            vals = parts[1:]
            if key not in meta:
                meta[key] = []
                order.append(key)
            meta[key] = vals
            n_samples = max(n_samples, len(vals))
    # pad
    for k in meta:
        if len(meta[k]) < n_samples:
            meta[k] += [""] * (n_samples - len(meta[k]))
    return meta, order, n_samples

def to_csv(meta, order, n_samples, gse, outdir):
    out = os.path.join(outdir, f"{gse}_samples.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        # header: sample index + each field as a column
        cols = ["sample_idx", "geo_sample_id"] + order
        w.writerow(cols)
        for i in range(n_samples):
            row = [i]
            # geo sample id from sample_geo_accession if present
            sid = ""
            if "sample_geo_accession" in meta:
                row_for = meta["sample_geo_accession"]
                sid = row_for[i] if i < len(row_for) else ""
            row.append(sid)
            for k in order:
                v = meta[k][i] if i < len(meta[k]) else ""
                v = re.sub(r"\s+", " ", v).strip()
                row.append(v)
            w.writerow(row)
    return out

if __name__ == "__main__":
    gses = sys.argv[1:] or ["GSE51618", "GSE209903", "GSE207937"]
    for gse in gses:
        p = os.path.join(GDRIVE, f"{gse}_series_matrix.txt.gz")
        if not os.path.exists(p):
            print(f"[skip] {gse} not found at {p}")
            continue
        print(f"=== parsing {gse} ===")
        meta, order, n = parse(p)
        out = to_csv(meta, order, n, gse, OUTDIR)
        print(f"  samples={n}  fields={len(order)}")
        print(f"  sample_id_field={'sample_geo_accession' in meta}")
        print(f"  -> {out}")
