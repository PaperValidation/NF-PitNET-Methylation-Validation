"""Permutation null for the cross-cohort nearest-neighbor mixing fraction.

Observed overall mixing (k=15) reported in manuscript: 0.3893
(per-cohort: Neou 0.51, GSE207937 0.51, GSE283928 0.22).

Null model: shuffle cohort labels across the FIXED PCA manifold (geometry
preserved) -> under H0 (cohort membership is meaningless) the mixing fraction
should centre on the size-weighted random expectation (~0.66). A significantly
lower observed value indicates real cohort cohesion (batch structure).
"""
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

SRC = r"C:\Users\tidus\PycharmProjects\Neurosurgery-NF-PitNET\results\Neou\joint_pca_coords.csv"
OUT = r"C:\Users\tidus\PycharmProjects\Neurosurgery-NF-PitNET\results\Neou\nn_mixing_permutation.json"

K = 15
P = 5000
RNG = 42

df = pd.read_csv(SRC)
cohorts = df["cohort"].values
X = df[[f"PC{i}" for i in range(1, 21)]].values.astype(float)
n = len(df)

# observed mixing fraction
nn = NearestNeighbors(n_neighbors=K + 1).fit(X)
_, idx = nn.kneighbors(X)
idx = idx[:, 1:]  # drop self
obs_fracs = np.zeros(n)
for i in range(n):
    neigh_cohort = cohorts[idx[i]]
    obs_fracs[i] = np.mean(neigh_cohort != cohorts[i])
observed_overall = float(obs_fracs.mean())
# per-cohort observed
per_cohort = {c: float(np.mean(obs_fracs[cohorts == c])) for c in np.unique(cohorts)}
# size-weighted random expectation
sizes = np.array([np.sum(cohorts == c) for c in np.unique(cohorts)])
random_expect = float((1 - sizes / n).dot(sizes / n))  # = mean over samples of P(diff cohort)

# permutation null
rng = np.random.default_rng(RNG)
perm_vals = np.empty(P)
for p in range(P):
    perm_cohort = rng.permutation(cohorts)
    f = np.zeros(n)
    for i in range(n):
        f[i] = np.mean(perm_cohort[idx[i]] != perm_cohort[i])
    perm_vals[p] = f.mean()

# one-sided lower tail: how extreme is observed (more separation) vs null
p_lower = (1 + np.sum(perm_vals <= observed_overall)) / (P + 1)
# two-sided
p_two = (1 + np.sum(np.abs(perm_vals - random_expect) >= np.abs(observed_overall - random_expect))) / (P + 1)

summary = {
    "k": K,
    "n_samples": n,
    "cohort_sizes": {c: int(np.sum(cohorts == c)) for c in np.unique(cohorts)},
    "observed_overall_mixing": round(observed_overall, 4),
    "per_cohort_observed_mixing": {c: round(v, 4) for c, v in per_cohort.items()},
    "size_weighted_random_expectation": round(random_expect, 4),
    "perm_n": P,
    "perm_null_mean": round(float(perm_vals.mean()), 4),
    "perm_null_sd": round(float(perm_vals.std()), 4),
    "perm_null_2p5": round(float(np.percentile(perm_vals, 2.5)), 4),
    "perm_null_97p5": round(float(np.percentile(perm_vals, 97.5)), 4),
    "p_value_one_sided_lower": p_lower,
    "p_value_two_sided": p_two,
}
import json
with open(OUT, "w") as f:
    json.dump(summary, f, indent=2)
print(json.dumps(summary, indent=2))
