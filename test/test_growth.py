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

    # Test which biomass components are producible by the model
    # It doesn't really test anyything, as in there is no way to pass or fail,
    # but it generates a heatmap of the producibility of the biomass components
    def test_biomass_component_producibility(self):
        # Load the TSV of the growth phenotypes
        growth_phenotypes = pd.read_csv(
            os.path.join(TESTFILE_DIR, "known_growth_phenotypes.tsv"), sep="\t"
        )

        # Load the model
        model = cobra.io.read_sbml_model("model.xml")

        # Filter the growth phenotypes to only include the carbon sources that it can grow on
        growth_phenotypes = growth_phenotypes[growth_phenotypes["growth"] == "Yes"]

        # Load the media definitions
        with open(
            os.path.join(TESTFILE_DIR, "media", "media_definitions.pkl"), "rb"
        ) as f:
            media_definitions = pickle.load(f)

        # Run the biomass producibility function on each of the models
        test_model(
            model, growth_phenotypes, media_definitions, biomass_rxn="bio1_biomass"
        )


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


# Broke out the function for actually testing the component producibility on a given medium
def test_medium(medium_dict, biomass_compounds, model):
    """Runs FBA for each biomass_compound as objective under `medium_dict`.
    Returns a dict {compound_id: True/False}."""
    results = {}
    # Set the medium
    model.medium = clean_media(model, medium_dict)

    # Loop over each biomass compound
    for cpd_id in biomass_compounds:
        # Get the human-friendly name of the metabolite
        cpd_name = model.metabolites.get_by_id(cpd_id).name

        # Set the objective to the sink/demand for cpd_id
        # e.g., "SK_cpd_id"
        model.objective = {model.reactions.get_by_id("SK_" + cpd_id): 1}

        # Optimize
        sol = model.optimize()

        # Check feasibility
        if (sol.status == "optimal") and (sol.objective_value > 1e-6):
            results[cpd_name] = True
        else:
            results[cpd_name] = False
    return results


def test_model(model, growth_phenotypes, media_definitions, biomass_rxn="bio1_biomass"):
    """
    Add exchange reactions for all metabolites in the model and test the producibility of the biomass components on the given media

    Args:
    model (cobra.Model): The model to test
    biomass_rxn (str): The ID of the biomass reaction in the model. Default is "bio1_biomass" (which is used in KBase models).
    """
    # Check that the growth phenotypes dataframe has the expected columns
    expected_columns = ["minimal_media", "c_source", "met_id", "growth"]
    if not all(col in growth_phenotypes.columns for col in expected_columns):
        raise ValueError(
            "growth_phenotypes dataframe must have columns "
            + ", ".join(expected_columns)
        )

    # Add sink reactions for all metabolites, but set the lower bound to 0
    # because by default the sink reactions are reversible, and so can be
    # used to import metabolites that are not in the media
    for metabolite in model.metabolites:
        # Check if there is already a sink reaction for this metabolite
        if "SK_" + metabolite.id not in [r.id for r in model.reactions]:
            model.add_boundary(metabolite, type="sink", lb=0)

    # Get the biomass composition from the model
    biomass_rxn = model.reactions.get_by_id(biomass_rxn)
    biomass_compounds = [
        met.id for met in biomass_rxn.metabolites if biomass_rxn.metabolites[met] < 0
    ]

    # Make a dictionary to store the producibility results
    biomass_producibility = {}

    # Add negative controls
    # TODO: Create the list of controls programmatically from the growth phenotypes data
    # Create some custom dicts for the "control" conditions:
    # 1) An empty medium
    empty_media = {}

    # 2) mbm minus carbon sources (assuming your minimal media has some "EX_CARBON" keys)
    #    This just filters out any exchange that might be for the primary carbon.
    mbm_no_carbon = {
        ex: lb
        for ex, lb in media_definitions["mbm_media"].items()
        if not ex.startswith("EX_") or "glc" not in ex.lower()
    }

    # 3) l1 minus carbon sources
    l1_no_carbon = {
        ex: lb
        for ex, lb in media_definitions["l1_media"].items()
        if not ex.startswith("EX_") or "glc" not in ex.lower()
    }

    # Combine the negative controls into a dictionary
    controls = {"empty": empty_media, "mbm_noC": mbm_no_carbon, "l1_noC": l1_no_carbon}

    # Test the negative controls
    for ctrl_name, ctrl_medium in controls.items():
        biomass_producibility[ctrl_name] = test_medium(
            ctrl_medium, biomass_compounds, model
        )

    # Check the producibility of the biomass componets on the different carbon sources
    for index, row in growth_phenotypes.iterrows():
        # Make an ID for the results that is combination of the minimal media name and the carbon source
        c_source = row["minimal_media"] + "_" + row["c_source"]
        # Make a dictionary to store the results for just this carbon source
        biomass_producibility[c_source] = {}
        # Set the model media to match the experimental media
        medium = media_definitions[row["minimal_media"]].copy()
        medium["EX_" + row["met_id"] + "_e0"] = (
            1000.0  # FIXME: I should set this to a consistent, lower value
        )
        # Test it
        biomass_producibility[c_source] = test_medium(medium, biomass_compounds, model)

    # Make a dataframe of the producibility results and save it to a CSV file
    df = pd.DataFrame.from_dict(biomass_producibility)
    # Save the dataframe to a CSV file and make the file name specific the the model.id
    df.to_csv(os.path.join(TESTFILE_DIR, model.id + "_biomass_producibility.csv"))

    # Plot the producibility results
    plot_prodcubility(model, df)


# Function for plotting the producibility results
def plot_prodcubility(model, df):
    # Convert the producibility results to a boolean
    df_binary = df.replace({True: 1, False: 0})

    # Define my own color palette
    my_palette = sns.color_palette(["red", "green"])

    # Plot the producibility results as a heatmap
    plt.figure(figsize=(10, 10))
    ax = sns.heatmap(
        df_binary,
        annot=False,
        cmap=my_palette,
        cbar_kws={"ticks": [0.25, 0.75], "label": "Producibile"},
    )

    # Now fix the colorbar tick labels
    cbar = ax.collections[0].colorbar
    cbar.set_ticklabels(["No", "Yes"])

    # Add white lines to separate the different cells
    for i in range(df_binary.shape[0]):
        plt.axhline(i, color="white", linewidth=0.5)
    for i in range(df_binary.shape[1]):
        plt.axvline(i, color="white", linewidth=1)

    # Make sure all y ticks/component names are shown
    plt.yticks(
        [x + 0.5 for x in range(df_binary.shape[0])], df_binary.index, fontsize=5
    )

    plt.xlabel("Biomass composition")
    plt.ylabel("Biomass compound")

    plt.tight_layout()

    # Save the plot
    plt.savefig(
        os.path.join(TESTFILE_DIR, model.id + "_biomass_producibility_heatmap.png")
    )


if __name__ == "__main__":
    unittest.main()
