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

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(SCRIPT_DIR))

# Load Prochlorococcus Genome-scale model
pro_cobra = cobra.io.read_sbml_model("iSO595v6.xml")
rename_pro_metabolites(pro_cobra)
model = c.model(pro_cobra)
model.initial_pop = [0, 0, 1e-7]
model.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

# Load the Alteromonas GEM and add it to the model as a second species
# amac_model = c.model(os.path.join(PROJECT_ROOT, "model.xml"))
# amac_model.initial_pop = [0, 0, 1e-7]
# amac_model.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

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

model.add_light("LightEX", absorption_biomass, absorption_water_680)

# Make layout with the COMETS toolbox
layout = c.layout(model)

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

# Run COMETS simulation
simulation = c.comets(layout, sim_params)
simulation.parameters.set_param("TotalBiomassLogName", "totalbiomasslog")
simulation.parameters.set_param("MediaLogName", "medialog")
simulation.run(delete_files=False)

# Remove COMETS temp files, keeping only the biomass/media logs and results
for pattern in [
    ".current_global_*",
    ".current_layout_*",
    ".current_script_*",
    ".current_package_*",
    "COBRAModel.cmd",
    "COMETS_manifest.txt",
]:
    for f in glob.glob(os.path.join(SCRIPT_DIR, pattern)):
        os.remove(f)

# Plot results
# Identify media metabolites with meaningful concentration changes
media = simulation.media.copy()
media["time"] = media["cycle"] * sim_params.all_params["timeStep"]
met_range = media.groupby("metabolite")["conc_mmol"].apply(lambda x: x.max() - x.min())
active_mets = met_range[met_range > 1e-10].index.tolist()

n_panels = 1 + len(active_mets)
fig, axes = plt.subplots(n_panels, 1, figsize=(10, 3 * n_panels), sharex=True)
if n_panels == 1:
    axes = [axes]

time_hours = simulation.total_biomass["cycle"] * sim_params.all_params["timeStep"]
axes[0].plot(time_hours, simulation.total_biomass.iloc[:, 1])
axes[0].set_ylabel("Biomass (gDW)")
axes[0].set_title("Prochlorococcus diel cycle simulation")

for ax, met in zip(axes[1:], active_mets):
    met_data = media[media["metabolite"] == met]
    ax.plot(met_data["time"], met_data["conc_mmol"])
    ax.set_ylabel("mmol")
    ax.set_title(f"{pro_cobra.metabolites.get_by_id(met).name} ({met})")

axes[-1].set_xlabel("Time (hours)")
plt.tight_layout()
plt.savefig(os.path.join(SCRIPT_DIR, "diel_cycle_results.png"), dpi=150)
