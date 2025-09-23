import os
import unittest

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch

# Set path to the `test_files` directory
TESTFILE_DIR = os.path.join(os.path.dirname(__file__), "test_files")

# Set path to the results directory
RESULTS_DIR = os.path.join(TESTFILE_DIR, "test_results")

# Define the colors that will go with presnt/absent
color_map = {"yes": "lightgreen", "no": "lightcoral"}


class TestPathways(unittest.TestCase):
    # Test that the reactions of specified pathways are present in the model.
    # Make a plot of the reaction presence/absence in the model.
    def test_pathway_reactions_present(self):
        # Load the model with cobrapy
        model = cobra.io.read_sbml_model("model.xml")

        # Load the CSV of the pathway reactions
        pathway_reactions = pd.read_csv(
            os.path.join(TESTFILE_DIR, "pathway_reactions.csv")
        )

        # Create a set of reaction IDs from the model for fast lookups
        model_reaction_ids = {r.id for r in model.reactions}

        # Use the vectorized `isin` method to get a boolean Series
        is_present = pathway_reactions["reaction_id"].isin(model_reaction_ids)

        # Map the boolean results to "Yes"/"No" strings
        pathway_reactions["in_model"] = is_present.map({True: "Yes", False: "No"})

        # Save the results to a new CSV file
        pathway_reactions.to_csv(
            os.path.join(RESULTS_DIR, "pathway_reaction_presence.csv"), index=False
        )

        # Plot a categorical heatmap of the reaction presence/absence for each pathway
        for pathway in pathway_reactions["pathway"].unique():
            plt_data = pathway_reactions[pathway_reactions["pathway"] == pathway].copy()
            plt_data = plt_data.pivot(
                index="reaction_id", columns="pathway", values="in_model"
            )
            plt.figure(figsize=(6, max(2, len(plt_data) * 0.3)))
            sns.heatmap(
                plt_data.replace({"Yes": 1, "No": 0}),
                annot=plt_data,
                cbar=False,
                vmin=0,
                vmax=1,
                cmap=[color_map["no"], color_map["yes"]],
                linewidths=0.5,
                linecolor="gray",
                fmt="",
            )
            plt.title(f"Reaction Presence in Model for {pathway} Pathway")
            plt.ylabel("Reaction ID")
            # Remove x-axis label and ticks
            plt.xlabel("")
            plt.xticks([], [])
            # Add a legend for the colors
            legend_elements = [
                Patch(facecolor=color_map["yes"], edgecolor="gray", label="Present"),
                Patch(facecolor=color_map["no"], edgecolor="gray", label="Absent"),
            ]
            plt.legend(
                handles=legend_elements,
                title="Reaction Presence",
                loc="center left",
                bbox_to_anchor=(1, 0.5),
            )
            plt.tight_layout()
            plt.savefig(
                os.path.join(RESULTS_DIR, f"pathway_reaction_presence_{pathway}.png")
            )
            plt.close()
