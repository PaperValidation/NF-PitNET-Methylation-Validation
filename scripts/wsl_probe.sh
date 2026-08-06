#!/usr/bin/env bash
# run inside WSL: wsl.exe -d Ubuntu-24.04 -- bash <thisfile>
B="/mnt/g/CMCC Sync/Medical Data/Neurosurgery/NF-PitNET"
echo "BASE=$B"
echo "=== downloaded files (size bytes, name) ==="
ls -la "$B" 2>/dev/null | awk 'NR>1 {print $5"\t"$9}'
echo "=== read test: line count GSE51618 matrix ==="
zcat "$B/GSE51618_series_matrix.txt.gz" 2>/dev/null | wc -l
echo "=== GSE283928 first header cols ==="
gunzip -c "$B/GSE283928_MatrixProcessed_Normalized.txt.gz" 2>/dev/null | head -1 | cut -f1-4
