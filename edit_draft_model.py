import os

import cobra
# Step 1: Make a draft model in KBase and download it

# Step 2 (This script):
# a. Add PHB to the biomass reaction
# b. Add all of Michelle's reactions to the model

REPO_DIR = os.path.dirname(os.path.abspath(__file__))
TESTFILE_DIR = os.path.join(REPO_DIR, "test", "test_files")

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