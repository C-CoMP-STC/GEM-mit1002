"""Make a scatter plot of the growth rate/yield/CUE for each single substrate
vs the abundance of that substrate in the cocktail."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.f2py import main
import pandas as pd
from scipy.stats import spearmanr
import seaborn as sns

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)


def main():
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

    # Merge the FBA results with their carbon concentrations
    singles = singles[singles["condition"] == "single"].copy()
    singles = singles.merge(
        top_10_exometabolites[["metabolite", "carbon_concentration"]],
        left_on="substrate",
        right_on="metabolite",
    )

    # Make a plot for the growth rate
    plot_y_vs_abundance(
        data=singles,
        col_name="growth_rate",
        y_label="Growth Rate on Single Substrate (1/hr)",
        title="Growth Rate vs Carbon Concentration in Cocktail",
        out_file_name="growth_rate_vs_carbon_concentration.png",
    )


def plot_y_vs_abundance(data, col_name, y_label, title, out_file_name):
    """Make a scatter plot of the growth rate/yield/CUE for each single substrate
    vs the abundance of that substrate in the cocktail."""
    # Calculate the Spearman correlation
    rho, pval = spearmanr(data["carbon_concentration"], data[col_name])

    # Make the figure
    fig, ax = plt.subplots(figsize=(8, 6))
    # Draw the scatterplot with different colors for each point/metabolite
    sns.scatterplot(
        data=data,
        x="carbon_concentration",
        y=col_name,
        hue="substrate",
        palette="tab10",
        s=100,
        edgecolor="black",
        linewidth=0.5,
    )
    # Draw a single regression line for all data (disable scatter to avaoid overlapping the colored points)
    sns.regplot(
        data=data,
        x="carbon_concentration",
        y=col_name,
        scatter=False,  # Disable the scatter points for the regression line
        line_kws={"color": "gray", "linestyle": "--"},
        ci=None,  # Don't plot the confidence interval for the regression line
    )
    # Annotate the plot with the Spearman correlation coefficient and p-value
    plt.text(
        0.05,
        0.95,
        f"Spearman ρ = {rho:.2f} (p = {pval:.2f})",
        transform=ax.transAxes,
        color="gray",
    )
    # Style
    plt.xlabel("Carbon Concentration in Cocktail (mM)", color="gray")
    plt.ylabel(y_label, color="gray")
    plt.title(title, color="gray")
    # Save
    plt.savefig(OUT_PATH / out_file_name, dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
