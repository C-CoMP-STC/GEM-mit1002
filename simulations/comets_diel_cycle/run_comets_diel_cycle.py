#!/Users/helenscott/opt/miniconda3/envs/med4-hot1a3/bin/python
import math
import os

import cometspy as c
import matplotlib.pyplot as plt

os.environ["COMETS_HOME"] = "/Applications/COMETS"

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

# Load Prochlorococcus Genome-scale model
model = c.model(os.path.join(SCRIPT_DIR, "iSO595v6.xml"))
model.initial_pop = [0, 0, 1e-7]
model.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

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
    "Ammonia[e]",
    "HCO3[e]",
    "CO2[e]",
    "H[e]",
    "Orthophosphate[e]",
    "H2O[e]",
    "Cadmium[e]",
    "Calcium_cation[e]",
    "Chloride_ion[e]",
    "Cobalt_ion[e]",
    "Copper[e]",
    "Fe2[e]",
    "Magnesium_cation[e]",
    "Molybdenum[e]",
    "K[e]",
    "Selenate[e]",
    "Sodium_cation[e]",
    "Strontium_cation[e]",
    "Sulfate[e]",
    "Zn2[e]",
    "Hydrogen_sulfide[e]",
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

# Run COMETS simulation
simulation = c.comets(layout, sim_params)
simulation.run(delete_files=False)

print(simulation.run_output)

# Plot results
fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

time_hours = simulation.total_biomass["cycle"] * sim_params.all_params["timeStep"]
axes[0].plot(time_hours, simulation.total_biomass.iloc[:, 1])
axes[0].set_ylabel("Biomass (gDW)")
axes[0].set_title("Prochlorococcus diel cycle simulation")

photon_media = simulation.media[simulation.media["metabolite"] == "Photon[e]"]
if not photon_media.empty:
    photon_time = photon_media["cycle"] * sim_params.all_params["timeStep"]
    axes[1].plot(photon_time, photon_media["conc_mmol"])
    axes[1].set_ylabel("Photon concentration (mmol)")
axes[1].set_xlabel("Time (hours)")

plt.tight_layout()
plt.savefig(os.path.join(SCRIPT_DIR, "diel_cycle_results.png"), dpi=150)
plt.show()
