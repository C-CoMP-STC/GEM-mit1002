"""Make the 9x9 pairwise epistasis heatmap."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)
pairs = pd.read_csv(OUT_PATH / "pairwise.csv")
singles = pd.read_csv(OUT_PATH / "single_substrate_results_total_c_60.csv")
singles = singles[singles["condition"] == "single"]
singles_lookup = dict(zip(singles["substrate"], singles["growth_rate"]))

# For each pair, compute the epistasis score:
# epistasis = growth(A+B) - (growth(A) + growth(B)) / 2
# (The /2 is because in the pair, each substrate is at half its single dose;
#  alternative is to compare to max or to a specific null model)

pairs["expected"] = pairs.apply(
    lambda r: (singles_lookup[r["substrate_a"]] + singles_lookup[r["substrate_b"]]) / 2,
    axis=1,
)
pairs["epistasis"] = pairs["growth_rate"] - pairs["expected"]

# Reshape into 9x9 matrix
substrates = sorted(set(pairs["substrate_a"]) | set(pairs["substrate_b"]))
matrix = pd.DataFrame(index=substrates, columns=substrates, dtype=float)
for _, row in pairs.iterrows():
    matrix.loc[row["substrate_a"], row["substrate_b"]] = row["epistasis"]
    matrix.loc[row["substrate_b"], row["substrate_a"]] = row["epistasis"]

# Plot with diverging colormap centered at zero
# (positive = synergy, negative = antagonism)
fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(
    matrix.values,
    cmap="RdBu_r",
    vmin=-matrix.abs().max().max(),
    vmax=matrix.abs().max().max(),
)
ax.set_xticks(range(len(substrates)))
ax.set_yticks(range(len(substrates)))
ax.set_xticklabels(substrates, rotation=45, ha="right")
ax.set_yticklabels(substrates)
fig.colorbar(im, label="Growth rate excess over expected")
plt.tight_layout()
fig.savefig(OUT_PATH / "figure_epistasis.png", dpi=150)
