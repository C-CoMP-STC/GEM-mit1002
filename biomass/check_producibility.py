import os
import pickle
import warnings

import cobra
import pandas as pd

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(FILE_DIR)
TESTFILE_DIR = os.path.join(REPO_DIR, "test", "test_files")


# Helper function for setting the media regardless if the exchange reaction is
# present in the model
# TODO: Move this to a helper file, and remove from the test growth file too
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
        # Set the objective to the sink/demand for cpd_id
        # e.g., "SK_cpd_id"
        model.objective = {model.reactions.get_by_id("SK_" + cpd_id): 1}

        # Optimize
        sol = model.optimize()

        # Check feasibility
        if (sol.status == "optimal") and (sol.objective_value > 1e-6):
            results[cpd_id] = True
        else:
            results[cpd_id] = False
    return results


# Load the model
model = cobra.io.read_sbml_model(
    os.path.join(REPO_DIR, "2025-01-08_Scott_draft-model-from-KBase.xml")
)

# Add exchange reactions for all metabolites
for metabolite in model.metabolites:
    model.add_boundary(metabolite, type="sink", lb=0)

# Load the growth pheonotype results
growth_phenotypes = pd.read_csv(
    os.path.join(TESTFILE_DIR, "known_growth_phenotypes.tsv"),
    sep="\t",
    header=0,
)

# Filter the growth phenotypes to only include the carbon sources that it can grow on
growth_phenotypes = growth_phenotypes[growth_phenotypes["growth"] == "Yes"]

# Load the media definitions
with open(os.path.join(TESTFILE_DIR, "media", "media_definitions.pkl"), "rb") as f:
    media_definitions = pickle.load(f)
minimal_media = media_definitions["minimal_media"]
mbm_media = media_definitions["mbm_media"]
l1_media = media_definitions["l1_media"]

# Get the biomass composition from the model
biomass_compounds = [
    met.id
    for met in model.reactions.bio1_biomass.metabolites
    if model.reactions.bio1_biomass.metabolites[met] < 0
]

# Check the producibility of the biomass componets on the different carbon sources
biomass_producibility = {}
for index, row in growth_phenotypes.iterrows():
    # Make an ID for the results that is combination of the minimal media name and the carbon source
    c_source = row["minimal_media"] + "_" + row["c_source"]
    # Make a dictionary to store the results for just this carbon source
    biomass_producibility[c_source] = {}
    # Set the model media to match the experimental media
    medium = eval(row["minimal_media"]).copy()
    medium["EX_" + row["met_id"] + "_e0"] = (
        1000.0  # FIXME: I should set this to a consistent, lower value
    )
    # Test it
    biomass_producibility[c_source] = test_medium(medium, biomass_compounds, model)

# Add negative controls
# Create some custom dicts for the "control" conditions:
# 1) An empty medium
empty_media = {}

# 2) mbm minus carbon sources (assuming your minimal media has some "EX_CARBON" keys)
#    This just filters out any exchange that might be for the primary carbon.
mbm_no_carbon = {
    ex: lb
    for ex, lb in mbm_media.items()
    if not ex.startswith("EX_") or "glc" not in ex.lower()
}

# 3) l1 minus carbon sources
l1_no_carbon = {
    ex: lb
    for ex, lb in l1_media.items()
    if not ex.startswith("EX_") or "glc" not in ex.lower()
}

# Combine the negative controls into a dictionary
controls = {"empty": empty_media, "mbm_noC": mbm_no_carbon, "l1_noC": l1_no_carbon}

# Test the negative controls
for ctrl_name, ctrl_medium in controls.items():
    biomass_producibility[ctrl_name] = test_medium(
        ctrl_medium, biomass_compounds, model
    )

# Make a dataframe of the producibility results and save it to a TSV file
df = pd.DataFrame.from_dict(biomass_producibility)
df.to_csv(os.path.join(FILE_DIR, "biomass_producibility.csv"))
