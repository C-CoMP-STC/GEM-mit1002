import os
import pickle
import sys

import cobra
import pandas as pd
from gem2cue import (
    utils,  # Import the working version (works with the med4-hot1a3 conda env)
)

# Set the output directory
OUT_DIR = os.path.dirname(os.path.realpath(__file__))

# Load the model and get the exchange reactions
model = cobra.io.read_sbml_model("model.xml")
c_ex_rxns = utils.get_c_ex_rxns(model)

# Load the results
with open(os.path.join(OUT_DIR, "results.pkl"), "rb") as f:
    cobra_results = pickle.load(f)

# Extract the carbon fate results from the FBA results and save them in a DataFrame
results_list = []
for key, fba_result in cobra_results.items():
    # Extract the carbon fates for the solution (both normalized and not normalized)
    c_fates = utils.extract_c_fates_from_solution(fba_result, c_ex_rxns, co2_ex_rxn='EX_cpd00011_e0', norm=False)
    uptake = c_fates[0]
    co2 = c_fates[1]
    organic_c = c_fates[2]
    biomass = c_fates[3]

    c_fates_norm = utils.extract_c_fates_from_solution(fba_result, c_ex_rxns, co2_ex_rxn='EX_cpd00011_e0', norm=True)
    co2_norm = c_fates_norm[0]
    organic_c_norm = c_fates_norm[1]
    biomass_norm = c_fates_norm[2]

    # Calculate CUE from the c fates (not using my function)
    cue = 1 - co2 / uptake

    # Calculate GGE from the c fates (not using my function)
    gge = 1 - (co2 + organic_c) / uptake

    results_list.append(
        {
            "sim_name": key,
            "oxygen_flux": fba_result.fluxes["EX_cpd00007_e0"],
            "uptake": uptake,
            "uptake_norm": 1,  # This is always 1 because the uptake is the reference
            "co2": co2,
            "co2_norm": co2_norm,
            "organic_c": organic_c,
            "organic_c_norm": organic_c_norm,
            "biomass": biomass,
            "biomass_norm": biomass_norm,
            "cue": cue,
            "gge": gge,
        }
    )

# Convert the results to a DataFrame
results = pd.DataFrame(results_list)

# Save the results
results.to_csv(os.path.join(OUT_DIR, "results.csv"), index=False)
