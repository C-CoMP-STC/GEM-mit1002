import os
import pickle
import sys
import warnings

import cobra
import pandas as pd
import matplotlib.pyplot as plt
from gem2cue import utils

# Import the plot styles (has global variables for colors)
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)
from plot_styles import *

# Set the output directory
OUT_DIR = os.path.dirname(os.path.realpath(__file__))

# Set a folder for the plots
output_folder = os.path.join(OUT_DIR, "plots")
# Check if the folder exists, if not, create it
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# ============================================================================
# SECTION 1: Load model and define media conditions
# ============================================================================

# Load in the ALT model using COBRApy
alt_cobra = cobra.io.read_sbml_model("model.xml")

# Make a medium with just glucose
# TODO: Use an uptake rate based on the NMR data
glc_medium_inf_o2 = {
    "EX_cpd00027_e0": 10,  # D-Glucose_e0
    "EX_cpd00007_e0": 1000,  # O2_e0
    # Remaining minimal media components
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# Make a medium with just acetate
# TODO: Use an uptake rate based on the NMR data
ace_medium_inf_o2 = {
    "EX_cpd00029_e0": 30,  # Acetate_e0
    "EX_cpd00007_e0": 1000,  # O2_e0
    # Remaining minimal media components
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# Make a medium with 2/3 glucose and 1/3 acetate
# FIXME: Need the equivalent amount of carbon available in the medium
glc_heavy_mix_medium_inf_o2 = {
    "EX_cpd00027_e0": 6.667,  # D-Glucose_e0
    "EX_cpd00029_e0": 10,  # Acetate_e0
    "EX_cpd00007_e0": 1000,  # O2_e0
    # Remaining minimal media components
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# Make a medium with 1/3 glucose and 2/3 acetate
# FIXME: Need the equivalent amount of carbon available in the medium
ace_heavy_mix_medium_inf_o2 = {
    "EX_cpd00027_e0": 3.333,  # D-Glucose_e0
    "EX_cpd00029_e0": 20,  # Acetate_e0
    "EX_cpd00007_e0": 1000,  # O2_e0
    # Remaining minimal media components
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# Make a medium with just glucose
# TODO: Use an uptake rate based on the NMR data
glc_medium_real_o2 = {
    "EX_cpd00027_e0": 10,  # D-Glucose_e0
    "EX_cpd00007_e0": 20,  # O2_e0
    # Remaining minimal media components
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# Make a medium with just acetate
# TODO: Use an uptake rate based on the NMR data
ace_medium_real_o2 = {
    "EX_cpd00029_e0": 30,  # Acetate_e0
    "EX_cpd00007_e0": 20,  # O2_e0
    # Remaining minimal media components
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# Make a medium with 2/3 glucose and 1/3 acetate
# FIXME: Need the equivalent amount of carbon available in the medium
glc_heavy_mix_medium_real_o2 = {
    "EX_cpd00027_e0": 6.667,  # D-Glucose_e0
    "EX_cpd00029_e0": 10,  # Acetate_e0
    "EX_cpd00007_e0": 20,  # O2_e0
    # Remaining minimal media components
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# Make a medium with 1/3 glucose and 2/3 acetate
# FIXME: Need the equivalent amount of carbon available in the medium
ace_heavy_mix_medium_real_o2 = {
    "EX_cpd00027_e0": 3.333,  # D-Glucose_e0
    "EX_cpd00029_e0": 20,  # Acetate_e0
    "EX_cpd00007_e0": 20,  # O2_e0
    # Remaining minimal media components
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

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


# Make a dictionary of the media with names
media = {
    "glc_medium_inf_o2": glc_medium_inf_o2,
    "ace_medium_inf_o2": ace_medium_inf_o2,
    "glc_heavy_mix_medium_inf_o2": glc_heavy_mix_medium_inf_o2,
    "ace_heavy_mix_medium_inf_o2": ace_heavy_mix_medium_inf_o2,
    "glc_medium_real_o2": glc_medium_real_o2,
    "ace_medium_real_o2": ace_medium_real_o2,
    "glc_heavy_mix_medium_real_o2": glc_heavy_mix_medium_real_o2,
    "ace_heavy_mix_medium_real_o2": ace_heavy_mix_medium_real_o2,
}

# ============================================================================
# SECTION 2: Run FBA and pFBA simulations
# ============================================================================

print("Running FBA and pFBA simulations...")
cobra_results = {}
for name, medium in media.items():
    # Run FBA
    alt_cobra.medium = clean_media(alt_cobra, medium)
    fba_result = alt_cobra.optimize()
    cobra_results[name + "_fba"] = fba_result
    # Run pFBA
    pfba_result = cobra.flux_analysis.pfba(alt_cobra)
    cobra_results[name + "_pfba"] = pfba_result

# ============================================================================
# SECTION 3: Calculate CUE and other metrics
# ============================================================================

print("Calculating CUE and related metrics...")

# Load the model and get the exchange reactions
c_ex_rxns = utils.get_c_ex_rxns(alt_cobra)

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
print(f"Results saved to {os.path.join(OUT_DIR, 'results.csv')}")

# ============================================================================
# SECTION 4: Generate plots
# ============================================================================

print("Generating plots...")

# Stacked bar plot of the carbon fates for the different conditions
data = results.set_index("sim_name")[["biomass", "organic_c", "co2"]]
g = carbon_fates_bar(data)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "carbon_fates.png"))
print(f"Saved plot: {os.path.join(output_folder, 'carbon_fates.png')}")

# Stacked bar plot of the normalized carbon fates for the different conditions
data_norm = results.set_index("sim_name")[["biomass_norm", "organic_c_norm", "co2_norm"]]
# Rename the columns to match the function
data_norm.columns = ["biomass", "organic_c", "co2"]
g_norm = carbon_fates_bar(data_norm)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "carbon_fates_norm.png"))
print(f"Saved plot: {os.path.join(output_folder, 'carbon_fates_norm.png')}")

# Subset the results to only include the pFBA on realistic O2
clean_data = data_norm[data_norm.index.str.contains("real_o2_pfba")]
g_clean = carbon_fates_bar(clean_data)
# Relabel the x tick labels
g_clean.set_xticklabels(
    ["Glucose", "Acetate", "Heavy Glucose Mix", "Heavy Acetate Mix"]
)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "carbon_fates_norm_clean.png"))
print(f"Saved plot: {os.path.join(output_folder, 'carbon_fates_norm_clean.png')}")

print("\nAnalysis complete!")