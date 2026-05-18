"""Make a scatter plot of the growth rate/yield/CUE for each single substrate
vs the abundance of that substrate in the cocktail."""

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


# Load the results from the simulations from the single + cocktail simulations
singles = pd.read_csv(
    TOP_10_DIR
    / "single_and_cocktail_sims"
    / "results"
    / "single_and_cocktail_results.csv"
)

# Plot the growth rate vs the carbon concentration for the single substrate condition
singles = singles[singles["condition"] == "single"].copy()
singles = singles.merge(
    top_10_exometabolites[["metabolite", "carbon_concentration"]],
    left_on="substrate",
    right_on="metabolite",
)
plt.figure(figsize=(6, 4))
sns.scatterplot(
    data=singles,
    x="carbon_concentration",
    y="growth_rate",
    color="gray",
    edgecolor="black",
    s=100,
)
plt.xlabel("Carbon Concentration in Cocktail (mM)", color="gray")
plt.ylabel("Growth Rate on Single Substrate (1/hr)", color="gray")
plt.title("Growth Rate vs Carbon Concentration of Substrate", color="gray")
plt.savefig(
    OUT_PATH / "growth_rate_vs_carbon_concentration.png", dpi=150, bbox_inches="tight"
)
