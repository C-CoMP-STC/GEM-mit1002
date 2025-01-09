import os

import cobra
import pandas as pd

# Load the model
model = cobra.io.read_sbml_model("2025-01-08_Scott_draft-model-from-KBase.xml")

# Add exchange reactions for all metabolites
for metabolite in model.metabolites:
    model.add_boundary(metabolite, type="sink")

###############################
# Get the biomass compositions
###############################
model_seed_path = "/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase/Templates"

# Load the biomass compositions from the TSV files
biomass_compositions = {}
biomass_compositions["core"] = pd.read_csv(os.path.join(model_seed_path, "Core", "BiomassCompounds.tsv"), sep="\t", header=0)['id'].tolist()
biomass_compositions["gram_positive"] = pd.read_csv(os.path.join(model_seed_path, "GramPositive", "BiomassCompounds.tsv"), sep="\t", header=0)['id'].tolist()
biomass_compositions["gram_negative"] = pd.read_csv(os.path.join(model_seed_path, "GramNegative", "BiomassCompounds.tsv"), sep="\t", header=0)['id'].tolist()

# Combine the biomass compositions, and keep only the unique values
biomass_compositions_all = set(biomass_compositions["core"] + biomass_compositions["gram_positive"] + biomass_compositions["gram_negative"])

# Check the producibility of the biomass compositions
biomass_producibility = {}
for biomass_compound in biomass_compositions_all:
    cpd_id = biomass_compound + "_c0"
    # Check if the compound is in the model
    try:
        model.metabolites.get_by_id(cpd_id)
    except KeyError:
        biomass_producibility[biomass_compound] = False
    # Check if the compound is producible
    else:
        # Set the objective to be the sink reaction for the compound
        model.objective = {model.reactions.get_by_id("SK_" + cpd_id): 1}
        # Run FBA
        sol = model.optimize()
        # Check if solution is feasible
        if sol.status != "optimal":
            biomass_producibility[biomass_compound] = False
        # Check if the compound is producible
        if sol.objective_value > 1e-6:
            biomass_producibility[biomass_compound] = True
        else:
            biomass_producibility[biomass_compound] = False

# Split the producibility results into the different biomass compositions
biomass_producibility_split = {}
biomass_producibility_split["core"] = {k: v for k, v in biomass_producibility.items() if k in biomass_compositions["core"]}
biomass_producibility_split["gram_positive"] = {k: v for k, v in biomass_producibility.items() if k in biomass_compositions["gram_positive"]}
biomass_producibility_split["gram_negative"] = {k: v for k, v in biomass_producibility.items() if k in biomass_compositions["gram_negative"]}

# Make a dataframe of the producibility results and save it to a TSV file
df = pd.DataFrame.from_dict(biomass_producibility_split)
df.to_csv("biomass/biomass_producibility.csv")
