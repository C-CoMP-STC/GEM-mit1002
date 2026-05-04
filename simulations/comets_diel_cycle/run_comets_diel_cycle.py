#!/Users/helenscott/opt/miniconda3/envs/med4-hot1a3/bin/python
import glob
import math
import os

import cobra
import cometspy as c
import matplotlib.pyplot as plt
import numpy as np
from pro_met_id_mapping import rename_pro_metabolites

os.environ["COMETS_HOME"] = "/Applications/COMETS"

FILE_DIR = os.path.dirname(os.path.realpath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(FILE_DIR))
RESULTS_DIR = os.path.join(FILE_DIR, "results")

# Make results directory if it doesn't exist
os.makedirs(RESULTS_DIR, exist_ok=True)

# Load Prochlorococcus Genome-scale model
pro_cobra = cobra.io.read_sbml_model("iSO595v6.xml")
rename_pro_metabolites(pro_cobra)
pro_model = c.model(pro_cobra)
pro_model.id = "iSO595v6"
pro_model.initial_pop = [0, 0, 1e-7]
pro_model.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

# Load the Alteromonas GEM and add it to the model as a second species
amac_model = c.model(os.path.join(PROJECT_ROOT, "model.xml"))
amac_model.initial_pop = [0, 0, 1e-7]
amac_model.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

# The ratio of chlorophyll is extracted from the model biomass-function
ci_dvchla = 0.016  # gr/gDW (Partensky 1993 / Casey 2016)
ci_dvchlb = 0.0013  # gr/gDW (Partensky 1993 / Casey 2016)
absorption_dvchla_680 = 0.0184  # m^2 mg^-1 (Bricaud et al., 2004)
absorption_dvchlb_680 = 0.0018  # m^2 mg^-1 (Bricaud et al., 2004)
absorption_water_680 = 0.465  # m^-1 (Pope and Fry, 1997)
wavelength = 680  # nm

diameter = 0.6  # um (Morel et al., 1993)
n_dash = (
    13.77 * 1e-3
)  # imaginary part of refractive index at 675 nm (Stramski et al. 2001)
size_parameter_alpha = (
    diameter * 1e3 * math.pi / wavelength
)  # The ratio between the cell size and wavelength
rho_dash = 4 * size_parameter_alpha * n_dash
Q_a = (
    1
    + (2 * math.exp(-rho_dash) / rho_dash)
    + 2 * (math.exp(-rho_dash) - 1) / rho_dash**2
)
packaging_effect = 1.5 * Q_a / rho_dash

# Calculate the Prochlorococcus specific biomass absorption coefficient in units m2/ g DW biomass
absorption_biomass = packaging_effect * (
    ci_dvchla * 1e3 * absorption_dvchla_680 + ci_dvchlb * 1e3 * absorption_dvchlb_680
)

pro_model.add_light("LightEX", absorption_biomass, absorption_water_680)

# Make an empty layout
layout = c.layout()
# Add the models to the layout
layout.add_model(pro_model)
layout.add_model(amac_model)

# Define medium
metabs = [
    "cpd00013_e0",  # Ammonia
    "HCO3[e]",
    "cpd00011_e0",  # CO2
    "cpd00067_e0",  # H+
    "cpd00009_e0",  # Phosphate
    "cpd00001_e0",  # H2O
    "Cadmium[e]",
    "cpd00063_e0",  # Ca2+
    "cpd00099_e0",  # Cl-
    "cpd00149_e0",  # Co2+
    "cpd00058_e0",  # Cu2+
    "cpd10515_e0",  # Fe+2
    "cpd00254_e0",  # Mg2+
    "Molybdenum[e]",
    "cpd00205_e0",  # K+
    "Selenate[e]",
    "cpd00971_e0",  # Na+
    "cpd09695_e0",  # Strontium
    "cpd00048_e0",  # Sulfate
    "cpd00034_e0",  # Zn2+
    "cpd00239_e0",  # HS- / bisulfide
]

for i in metabs:
    layout.set_specific_metabolite(i, 1000)
    layout.set_specific_static(i, 1000)

# Set light conditions by defining parameters
layout.set_global_periodic_media(
    metabolite="Photon[e]",
    function="half_sin",
    amplitude=0.04,
    period=24,
    phase=0,
    offset=0,
)

# Set simulation parameters
sim_params = c.params()
sim_params.all_params["maxCycles"] = 480
sim_params.all_params["timeStep"] = 0.1
sim_params.all_params["defaultDiffConst"] = 0
sim_params.all_params["writeMediaLog"] = True
sim_params.all_params["MediaLogRate"] = 1
sim_params.all_params["writeFluxLog"] = True
sim_params.all_params["FluxLogRate"] = 1

# Run COMETS simulation
simulation = c.comets(layout, sim_params)
simulation.run()

# Save the simulation results to a CSV file
simulation.total_biomass.to_csv(
    os.path.join(RESULTS_DIR, "total_biomass.csv"), index=False
)
simulation.media.to_csv(os.path.join(RESULTS_DIR, "media_log.csv"), index=False)
for model_id, flux_df in simulation.fluxes_by_species.items():
    flux_df["time"] = flux_df["cycle"] * sim_params.all_params["timeStep"]
    flux_df.to_csv(os.path.join(RESULTS_DIR, f"fluxlog_{model_id}.csv"), index=False)
