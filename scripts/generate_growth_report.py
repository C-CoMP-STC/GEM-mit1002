import os

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from gem_utilities import biomass

# Define paths relative to the script or project root
# It's better practice to define a project root
PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))

import sys

sys.path.insert(0, str(PROJECT_ROOT))

from tools.media import MEDIA  # noqa: E402
from tools.phenotypes import evaluate_phenotypes  # noqa: E402

DATA_DIR = os.path.join(PROJECT_ROOT, "data")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "scripts", "results")

# Load the media definitions
media_definitions = MEDIA

# Define a dictionary of human-friendly versions of the media names
# The key is the name in the media_definitions dictionary
# The value is the human-friendly name to use the table
media_names = {
    "l1": "L1",
    "mbm": "Minimal Basal Medium (Moran Lab)",
    "promm_no_c": "ProMM",
    "marine_broth_wo_yeast_and_peptone": "Marine Broth",
    "marine_broth_wo_yeast_and_peptone_no_n": "Marine Broth (No Nitrogen)",
    "swm": "Seawater Medium",
}


def generate_growth_phenotype_report(model: cobra.Model):
    # Evaluate every condition with the shared scorer in tools/phenotypes.py,
    # the same one used by test/test_growth.py and by the figure 2 pipeline in
    # curation_process/run_tests_on_prs.py.
    #
    # This function used to carry its own copy of the simulation loop, which
    # disagreed with that scorer in three ways:
    #   * it added a row's metabolites to the medium only when *every* one of
    #     them had an exchange reaction, so "Methionine, Pyruvate" was simulated
    #     on a carbon-free medium rather than on the pyruvate;
    #   * it never checked the solution status, and the model's forced ATP
    #     maintenance bound makes a medium with no usable carbon infeasible
    #     rather than zero-growth, so those rows reported whatever the solver
    #     had left in the primal -- the negative "growth rates" in the table;
    #   * it therefore could not distinguish an infeasible solve from a genuine
    #     prediction of no growth.
    # See the module docstring of tools/phenotypes.py for the full rationale.
    growth_phenotypes = evaluate_phenotypes(model)

    # Re-expose the two columns under the names the heatmap below expects.
    # An unevaluable solve has predicted=None; render it as "Unsure" (gray)
    # rather than silently as no growth.
    growth_phenotypes["all_ex_rxn_present"] = [
        "No" if missing else "Yes"
        for missing in growth_phenotypes["missing_exchanges"]
    ]
    growth_phenotypes["pred_growth"] = growth_phenotypes["predicted"].fillna("Unsure")

    # Beautify and save the table
    beautify_table(growth_phenotypes)

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

    # Make a dictionary for the phenotypes to numbers
    value_to_int = {"Unsure": 0, "No": 1, "Yes": 2}
    n = len(value_to_int)

    # Create an annotation data frame for the text labels on the heatmap
    annotation_key = {"No": "No Exchange", "Yes": ""}
    annot_df = (
        growth_phenotypes["all_ex_rxn_present"].replace(annotation_key).to_frame()
    )
    annot_df.rename(columns={"all_ex_rxn_present": "FBA"}, inplace=True)
    annot_df["Experimental"] = ""
    # Sort the columns to match the order of the heatmap
    annot_df = annot_df[["Experimental", "FBA"]]

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

    # Make a colormap of specified colors (in numerical order for the phenotypes)
    # cmap = ['gray', '#F18F01', '#399E5A'] # Gray, orange, green
    cmap = ["#5E5E5E", "#FF7D0A", "#024064"]  # C-CoMP gray, orange, and dark blue

    # Dynamically set the figure height based on the number of rows
    fig_height = max(10, len(growth_phenotypes) * 0.4)  # 0.4 inches per row
    # Plot the heatmap
    # Use constrained_layout to prevent cutting off y-axis/colorbar labels
    fig, ax = plt.subplots(
        figsize=(8, fig_height),
        constrained_layout=True,
    )
    sns.heatmap(
        growth_phenotypes.replace(value_to_int),
        cmap=cmap,
        linewidths=4,
        linecolor="white",
        annot=annot_df,
        fmt="",
        annot_kws={"fontsize": 8},  # Smaller font size for annotation
        ax=ax,
    )

    # Modify colorbar:
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
        top=True,
        labeltop=True,
    )

    # Make sure that every y-tick is shown
    ax.set_yticks([i + 0.5 for i in range(len(growth_phenotypes))])
    ax.set_yticklabels(growth_phenotypes.index, rotation=0)

    # Save the figure
    plt.savefig(os.path.join(RESULTS_DIR, "exp_vs_pred_growth_phenotypes.png"))


def beautify_table(exp_pred_table: pd.DataFrame):
    # Get all of the unique minimal media in the table
    unique_minimal_media = exp_pred_table["minimal_media"].unique()
    # Make a dictionary where the keys are the minimal media and the values are the nitrogen-containing compounds in that media, separated by commas
    n_dict = {}
    # For each unique minimal media, look up the nitrogen-containing compounds in that media and add a column to the table with that information
    for minimal_media in unique_minimal_media:
        # Get the media definition for that minimal media
        media_def = media_definitions[minimal_media]
        # Get the nitrogen-containing compounds in that media definition
        n_containing_compounds = []
        for ex_rxn in media_def:
            # Get the metabolite ID from the reaction ID
            # (remove the "EX_" prefix)
            met_id = ex_rxn[3:]  # Remove "EX_"
            # Try to get the metabolite object from the model
            try:
                met = model.metabolites.get_by_id(met_id)
            except KeyError:
                continue
            # Check if the metabolite contains nitrogen in its formula
            if "N" in met.elements:
                # Ignore vitamins (thiamin and vitamin B12)
                if met.id in ["cpd00305_e0", "cpd03424_e0"]:
                    continue
                # If it does, add the name of the metabolite to the list of nitrogen-containing compounds
                # Remove the " [e0]" from the end of the name, if it exists
                if met.name.endswith(" [e0]"):
                    met_name = met.name[:-5]
                else:
                    met_name = met.name
                n_containing_compounds.append(met_name)
        # Store the list of nitrogen-containing compounds for this minimal media
        # as a string separated by commas
        n_dict[minimal_media] = ", ".join(n_containing_compounds)
    # Add a column to the table with the nitrogen-containing compounds, separated by commas
    exp_pred_table["Medium N Source(s)"] = exp_pred_table["minimal_media"].map(n_dict)

    # Replacing the media names in the minimal_media column with the human-friendly versions
    exp_pred_table["minimal_media"] = exp_pred_table["minimal_media"].map(media_names)

    # Round growth rates below (absolute value) 1e-3 to 0
    # To avoid negative 0s
    exp_pred_table["fba_growth_rate"] = exp_pred_table["fba_growth_rate"].apply(
        lambda x: 0 if abs(x) < 1e-3 else x
    )

    # Round the FBA predicted growth rates to 3 decimal places for easier readability
    exp_pred_table["fba_growth_rate"] = exp_pred_table["fba_growth_rate"].round(3)

    # Subset the columns we want
    exp_pred_table = exp_pred_table[
        [
            "minimal_media",
            "Medium N Source(s)",
            "c_source",
            "growth",
            "reference",
            "fba_growth_rate",
        ]
    ]
    # Sort the table by the minimal media and then the carbon source
    exp_pred_table = exp_pred_table.sort_values(by=["minimal_media", "c_source"])
    # And rename the columns to be more descriptive
    exp_pred_table = exp_pred_table.rename(
        columns={
            "minimal_media": "Minimal Media",
            "c_source": "Added Metabolite(s)",
            "growth": "Experimental Growth",
            "reference": "Reference",
            "fba_growth_rate": "FBA Predicted Growth Rate",
        }
    )
    # Save as a tsv
    exp_pred_table.to_csv(
        os.path.join(RESULTS_DIR, "known_growth_phenotypes_w_pred.tsv"),
        index=False,
        sep="\t",
    )
    # Save as excel
    exp_pred_table.to_excel(
        os.path.join(RESULTS_DIR, "known_growth_phenotypes_w_pred.xlsx"), index=False
    )


def generate_biomass_producibility_report(model: cobra.Model):
    # Load the TSV of the growth phenotypes
    growth_phenotypes = pd.read_csv(
        os.path.join(DATA_DIR, "known_growth_phenotypes.tsv"),
        sep="\t",
        converters={"met_id": lambda x: x.split(",")},
    )

    # Filter the growth phenotypes to only include the carbon sources that it can grow on
    growth_phenotypes = growth_phenotypes[growth_phenotypes["growth"] == "Yes"]

    # Run the biomass producibility function on each of the models
    sink_options = [False, True]
    for option in sink_options:
        biomass.check_biomass_producibility(
            model,
            growth_phenotypes,
            media_definitions,
            sinks_for_all=option,
            out_dir=RESULTS_DIR,
        )


if __name__ == "__main__":
    # Ensure the results directory exists
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(PROJECT_ROOT, "model.xml"))

    # Generate the reports
    generate_growth_phenotype_report(model)
    generate_biomass_producibility_report(model)
