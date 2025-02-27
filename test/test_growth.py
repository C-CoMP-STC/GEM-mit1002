import os
import pickle
import unittest
import warnings

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set path to the `test_files` directory
TESTFILE_DIR = os.path.join(os.path.dirname(__file__), "test_files")

# Load the media definitions
with open(os.path.join(TESTFILE_DIR, "media", "media_definitions.pkl"), "rb") as f:
    media_definitions = pickle.load(f)
minimal_media = media_definitions["minimal"]


class TestGrowthPhenotypes(unittest.TestCase):
    # Check that there is no growth on a media with no carbon sources
    def test_growth_w_0_C(self):
        # Load the model with cobrapy
        model = cobra.io.read_sbml_model("model.xml")

        # Set the media so that there are no carbon sources
        model.medium = clean_media(model, minimal_media)

        # Run the model
        sol = model.optimize()

        # Check that no biomass is produced
        self.assertEqual(sol.objective_value, 0)

    # Test that there is growth on the expected pheontypes, and that there
    # is no growth on the expected phenotypes. Make a plot of the growth
    # phenotypes.
    def test_expected_growth_phenotypes(self):
        # Load the TSV of the growth phenotypes
        growth_phenotypes = pd.read_csv(
            os.path.join(TESTFILE_DIR, "known_growth_phenotypes.tsv"), sep="\t"
        )

        # Load the model
        model = cobra.io.read_sbml_model("model.xml")

        # Loop through the growth phenotpes, and add the carbon source to the
        # minimal media, run FBA and check if the model grows
        pred_growth = []
        for index, row in growth_phenotypes.iterrows():
            minimal_media = media_definitions[row["minimal_media"]].copy()
            # Check if the model has an exchange reaction for the metabolite
            if "EX_" + row["met_id"] + "_e0" in [r.id for r in model.reactions]:
                # If it does, add the exchange reaction to the minimal media used
                minimal_media["EX_" + row["met_id"] + "_e0"] = 1000.0
            else:
                # If it doesn't have the exchange reaction, add "No Exchange"
                pred_growth.append("No Exchange")
                # Give a warning if growth was expected
                if row["growth"] == "Yes":
                    warnings.warn(
                        "Model did not have an exchange reaction for "
                        + row["c_source"]
                        + ", but growth was expected."
                    )
                continue
            # Set the media
            model.medium = clean_media(model, minimal_media)
            # Run the model
            sol = model.optimize()
            # Check if the model grows
            if sol.objective_value > 0:
                # If it does, add 'Y' to the list
                pred_growth.append("Yes")
                # Give a warning if growth was not expected
                if row["growth"] == "No":
                    warnings.warn(
                        "Model grew on "
                        + row["c_source"]
                        + ", but growth was not expected."
                    )
            else:
                # If it doesn't, add 'N' to the list
                pred_growth.append("No")
                # Give a warning if growth was expected
                if row["growth"] == "Yes":
                    warnings.warn(
                        "Model did not grow on "
                        + row["c_source"]
                        + ", but growth was expected."
                    )

        # Add the list as a new column in the dataframe
        growth_phenotypes["pred_growth"] = pred_growth

        # Save the dataframe as a TSV
        growth_phenotypes.to_csv(
            os.path.join(TESTFILE_DIR, "known_growth_phenotypes_w_pred.tsv"),
            sep="\t",
            index=False,
        )

        # Plot a categorical heatmap of the growth phenotypes, where the rows
        # are the metabolites and the columns are the experimental and predicted
        # growth phenotypes. Show growth as blue and no growth as orange, and
        # unsure as gray
        # First, make a new dataframe with the metabolites as the rows and the
        # experimental and predicted growth phenotypes as the columns
        # Combine the values of "minimal_media" and "c_source" into one column
        growth_phenotypes["c_source"] = (
            growth_phenotypes["minimal_media"] + " " + growth_phenotypes["c_source"]
        )
        # And set it as the index
        growth_phenotypes = growth_phenotypes.set_index("c_source")
        # Subset the other columns, to have just the growth and predicted growth
        growth_phenotypes = growth_phenotypes[["growth", "pred_growth"]]

        # Rename the columns and the index to be longer/more descriptive
        growth_phenotypes.index.name = "Media/Carbon Source"
        growth_phenotypes = growth_phenotypes.rename(
            columns={
                "growth": "Experimental",
                "pred_growth": "FBA",
            }
        )

        # Replace all of the "No Exchange" values with "No"
        growth_phenotypes = growth_phenotypes.replace("No Exchange", "No")

        # Make a dictionary for the phenotypes to numbers
        value_to_int = {"Unsure": 0, "No": 1, "Yes": 2}
        n = len(value_to_int)

        # Make a colormap of specified colors (in numerical order for the phenotypes)
        # cmap = ['gray', '#F18F01', '#399E5A'] # Gray, orange, green
        cmap = ["#5E5E5E", "#FF7D0A", "#024064"]  # C-CoMP gray, orange, and dark blue

        # Plot the heatmap
        ax = sns.heatmap(
            growth_phenotypes.replace(value_to_int),
            cmap=cmap,
            linewidths=4,
            linecolor="white",
        )

        # modify colorbar:
        colorbar = ax.collections[0].colorbar
        r = colorbar.vmax - colorbar.vmin
        colorbar.set_ticks([colorbar.vmin + r / n * (0.5 + i) for i in range(n)])
        colorbar.set_ticklabels(list(value_to_int.keys()))

        # Move the x-axis labels to the top
        plt.tick_params(
            axis="both",
            which="major",
            labelsize=10,
            labelbottom=False,
            bottom=False,
            top=False,
            labeltop=True,
        )

        # Make sure that every y-tick is shown
        ax.set_yticks([i + 0.5 for i in range(len(growth_phenotypes))])
        ax.set_yticklabels(growth_phenotypes.index, rotation=0)

        # Make sure that the y-axis labels are not cut off
        plt.tight_layout()

        # Save the figure
        plt.savefig(os.path.join(TESTFILE_DIR, "exp_vs_pred_growth_phenotypes.png"))


# Helper function for setting the media regardless if the exchange reaction is
# present in the model
# TODO: Move this to a helper file
def clean_media(model, media):
    """clean_media
    Removes exchange reactions from the media that are not present in the model

    Parameters
    ----------
    model : cobra.Model
        The model to set the media for.
    media : dict
        A dictionary where the keys are the exchange reactions for the metabolites
        in the media, and the values are the lower bound for the exchange reaction.

    Returns
    -------
    dict
        A dictionary where the keys are the exchange reactions for the metabolites
        in the media, and the values are the lower bound for the exchange reaction
    """
    # Make an empty dictionary for the media
    clean_medium = {}
    # Loop through the media and set the exchange reactions that are present
    for ex_rxn, lb in media.items():
        if ex_rxn in [r.id for r in model.reactions]:
            clean_medium[ex_rxn] = lb
        else:
            warnings.warn(
                "Model does not have the exchange reaction "
                + ex_rxn
                + ", so it was not set in the media."
            )

    # Return the clean medium
    return clean_medium


if __name__ == "__main__":
    unittest.main()
