#!/usr/bin/env python3
"""Quick NaN audit for Neou beta matrix."""
import numpy as np, os

from _paths import DATASET_PATH as BASE
NEU_BIN = os.path.join(BASE, "E-MTAB-7762", "beta_out", "neou_beta_full.bin")

with open(NEU_BIN, "rb") as f:
    nrow = int(np.fromfile(f, dtype=np.int32, count=1)[0])
    ncol = int(np.fromfile(f, dtype=np.int32, count=1)[0])
    data = np.fromfile(f, dtype=np.float64, count=nrow * ncol)
Bmat = data.reshape(nrow, ncol, order="F")  # probes x samples

total = Bmat.size
nan_count = np.isnan(Bmat).sum()
print(f"Neou beta matrix: {nrow} probes x {ncol} samples = {total} values")
print(f"NaN count: {nan_count} ({100*nan_count/total:.4f}%)")
print(f"NaN probes (all-NaN): {int(np.isnan(Bmat).all(axis=1).sum())}")
print(f"NaN samples (any NaN): {int(np.isnan(Bmat).any(axis=0).sum())}")

# per-sample NaN fraction
per_samp_nan = np.isnan(Bmat).mean(axis=0)
print(f"Per-sample NaN fraction: min={per_samp_nan.min():.6f} max={per_samp_nan.max():.6f} mean={per_samp_nan.mean():.6f}")
# any sample with >5% NaN?
high_nan = np.where(per_samp_nan > 0.05)[0]
print(f"Samples with >5% NaN: {len(high_nan)}")
if len(high_nan):
    print(f"  indices: {high_nan.tolist()}")
