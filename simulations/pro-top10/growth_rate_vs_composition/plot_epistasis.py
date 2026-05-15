"""Make the 9x9 pairwise epistasis heatmap."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)

# Load the original metabolite list to get their concentrations to sort by
TOP_10_PATH = TOP_10_DIR / "data" / "top10_exometabolites.csv"
top_10_exometabolites = pd.read_csv(TOP_10_PATH)
# Filter the data to only include the Prochlorococcus marinus metabolites
top_10_exometabolites = top_10_exometabolites[
    top_10_exometabolites["organism"] == "Prochlorococcus marinus"
].copy()
# Remove the metabolites that don't support growth/that I didn't test
# In other scripts I use the model to check this, but I don't want to have to
# load the model here just for that, and I know from before that it's just Met
# and Phe that don't have exchanges
top_10_exometabolites = top_10_exometabolites[
    ~top_10_exometabolites["metabolite"].isin(["methionine", "phenylalanine"])
].copy()
# Sort that list by the carbon concentration
top_10_exometabolites = top_10_exometabolites.sort_values(
    by="carbon_concentration", ascending=False
).copy()

# Get the desired order of substrates
substrate_order = top_10_exometabolites["metabolite"].tolist()

# Load the results from the pairwise simulations and the single substrate
# simulations from the single + cocktail simulations
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
matrix = pd.DataFrame(index=substrate_order, columns=substrate_order, dtype=float)
for _, row in pairs.iterrows():
    matrix.loc[row["substrate_a"], row["substrate_b"]] = row["epistasis"]
    matrix.loc[row["substrate_b"], row["substrate_a"]] = row["epistasis"]

# Create a mask for the upper triangle
mask = np.triu(np.ones_like(matrix, dtype=bool))

# Plot with diverging colormap centered at zero
# (positive = synergy, negative = antagonism)
fig, ax = plt.subplots(figsize=(8, 7))
sns.heatmap(
    matrix,
    mask=mask,  # Apply the mask
    cmap="RdBu_r",
    vmax=matrix.abs().max().max(),
    vmin=-matrix.abs().max().max(),
    annot=True,  # Optionally, show the values in the cells
    fmt=".2f",  # Format annotations to 2 decimal places
    linewidths=0.5,
    ax=ax,
)

ax.set_title("Pairwise Growth Epistasis")
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

fig.tight_layout()
fig.savefig(OUT_PATH / "figure_epistasis.png", dpi=150)
