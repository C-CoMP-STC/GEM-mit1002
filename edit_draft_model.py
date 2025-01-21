import ast
import json
import os

import cobra
import pandas as pd

# Step 1: Make a draft model in KBase and download it

# Step 2 (This script):
# a. Add PHB to the biomass reaction
# b. Add all of Michelle's reactions to the model

REPO_DIR = os.path.dirname(os.path.abspath(__file__))

# Load the draft model from KBase
model = cobra.io.read_sbml_model(
    os.path.join(REPO_DIR, "2025-01-08_Scott_draft-model-from-KBase.xml")
)

##################################################
# Add PHB to the biomass reaction
##################################################
# Add the metabolite to the model
phb = cobra.Metabolite(
    "cpd27893_c0",
    formula="C12H18O7R",
    name="Poly-Hydroxybutyrate",
    compartment="c0",
)
model.add_metabolites([phb])

# Find the biomass reaction
biomass_rxn = model.reactions.get_by_id("bio1_biomass")
# Add PHB to the biomass reaction
# FIXME: Does that mean that the biomass reaction is not balanced?
biomass_rxn.add_metabolites(
    {
        model.metabolites.get_by_id("cpd27893_c0"): -0.1,
    }
)

# Save the model
cobra.io.write_sbml_model(
    model,
    os.path.join(REPO_DIR, "2025-01-08_Scott_draft-model-from-KBase_PHB.xml"),
)

##################################################
# Add Michelle's reactions to the model
##################################################
# Load the database of Michelle's reactions with the ModelSEED IDs
michelle_rxns = pd.read_csv(
    os.path.join(REPO_DIR, "Pangenome from Michelle", "database_w_MNX_SEED.csv"),
    header=0,
)

# Filter the list to only include ones with a "1" for "Inferred Presence" and with something in "ModelSEED ID"
# Make an explicit copy after filtering to avoid SettingWithCopyWarning
rxns_to_add = michelle_rxns[
    (michelle_rxns["Inferred Presence"] == 1) & (michelle_rxns["ModelSEED ID"].notnull())
].copy()

# Check if any of the "ModelSEED ID" values are empty
# rxns_to_add["ModelSEED ID"].isnull().sum()

# Convert the "ModelSEED ID" column from a string to a list
rxns_to_add["ModelSEED ID"] = rxns_to_add["ModelSEED ID"].apply(ast.literal_eval)

# Combine all of the values in the "ModelSEED ID" column into a list
lists_of_rxn_ids = rxns_to_add["ModelSEED ID"].tolist()
rxn_ids = [
    x
    for item in lists_of_rxn_ids
    for x in (item if isinstance(item, list) else [item])
]
# Get just the unique values
rxn_ids = list(set(rxn_ids))

# Check if the ModelSEED IDs are obsolete
# Need to load in the ModelSEED database first
modelseed_db = json.load(
    open(
        "/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase/Biochemistry/reactions.json"
    )
)
# Convert to a dictionary with the ModelSEED IDs as the keys for easier searching
modelseed_db = {rxn["id"]: rxn for rxn in modelseed_db}

# Subset the rxn_ids list to only include non-obsolete reactions
rxn_ids_to_add = [rxn_id for rxn_id in rxn_ids if rxn_id in modelseed_db and not modelseed_db[rxn_id]["is_obsolete"]]

# Add the remaining, non-obsolete reactions to the model
# TODO: Annotate to say why these reactions are being added
for rxn_id in rxn_ids_to_add:
    rxn = modelseed_db[rxn_id]
    reaction = cobra.Reaction(
        rxn["id"] + '_c0',
        name=rxn["name"],
        # TODO: Handle the reversibility
    )
    # Add the metabolites to the reaction
    for met in rxn["stoichiometry"].split(";"):
        coeff = met.split(":")[0]
        met_id_no_cmpt = met.split(":")[1]
        if met.split(":")[2] == '0':
            met_id = met_id_no_cmpt + '_c0'
        elif met.split(":")[2] == '2':
            print("CHECK COMPARTMENT")
            met_id = met_id_no_cmpt + '_e0'
        if met_id not in model.metabolites:
            met_obj = cobra.Metabolite(
                met_id,
                name=met_id,
                compartment=met_id.split('_')[1],
            )
            model.add_metabolites([met_obj])
        met_obj = model.metabolites.get_by_id(met_id)
        reaction.add_metabolites({met_obj: float(coeff)})
    model.add_reactions([reaction])

# Save the model
cobra.io.write_sbml_model(
    model,
    os.path.join(REPO_DIR, "2025-01-08_Scott_draft-model-from-KBase_PHB_Michelle.xml"),
)

# Step 3: Upload the model to KBase and do gap filling
