import os

import cobra
import matplotlib.pyplot as plt
import pandas as pd

FILE_DIR = os.path.dirname(os.path.realpath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(FILE_DIR))

TIME_STEP = 0.1  # Make sure this matches the time step used in the COMETS simulation

# Load the simulation results
biomass = pd.read_csv(os.path.join(FILE_DIR, "total_biomass.csv"))
media = pd.read_csv(os.path.join(FILE_DIR, "media_log.csv"))

# Load the models (to get metabolite names)
pro_cobra = cobra.io.read_sbml_model("iSO595v6.xml")
amac_cobra = cobra.io.read_sbml_model(os.path.join(PROJECT_ROOT, "model.xml"))

# Get a lookup dinctionary of metabolite IDs and names across both models
met_names = {met.id: met.name for met in pro_cobra.metabolites}
met_names.update({met.id: met.name for met in amac_cobra.metabolites})

# Add "time" columns to each dataframe based on the "cycle" column and the time step, to make plotting easier
biomass["time"] = biomass["cycle"] * TIME_STEP
media["time"] = media["cycle"] * TIME_STEP

# Identify media metabolites with meaningful concentration changes
met_range = media.groupby("metabolite")["conc_mmol"].apply(lambda x: x.max() - x.min())
active_mets = met_range[met_range > 1e-10].index.tolist()

# Make one plot with a subplot for each active metabolite's concentration over time
n_panels = len(active_mets)
fig, axes = plt.subplots(n_panels, 1, figsize=(10, 3 * n_panels), sharex=True)
if n_panels == 1:
    axes = [axes]
for ax, met in zip(axes[0:], active_mets):
    met_data = media[media["metabolite"] == met]
    ax.plot(met_data["time"], met_data["conc_mmol"])
    ax.set_ylabel("mmol")
    ax.set_title(f"{met_names.get(met, met)} ({met})")
axes[-1].set_xlabel("Time (hours)")
plt.tight_layout()
plt.savefig(os.path.join(FILE_DIR, "diel_cycle_individual_mets.png"), dpi=150)

# Make a plot with biomass (for both models) and light on twin y-axes
fig, ax1 = plt.subplots(figsize=(10, 6))
ax1.plot(biomass["time"], biomass["COBRAModel"], label="Prochlorococcus", color="green")
ax1.plot(biomass["time"], biomass["iHS4156"], label="Alteromonas", color="blue")
ax1.set_xlabel("Time (hours)")
ax1.set_ylabel("Biomass (gDW)")
ax1.legend(loc="upper left")
ax2 = ax1.twinx()
photon = media[media["metabolite"] == "Photon[e]"]
ax2.plot(
    photon["time"],
    photon["conc_mmol"],
    label="Light",
    color="orange",
    linestyle="--",
)
ax2.set_ylabel("Light intensity (umol photons m^-2 s^-1)")
ax2.legend(loc="upper right")
plt.title("Biomass and Light over Diel Cycle")
plt.tight_layout()
plt.savefig(os.path.join(FILE_DIR, "diel_cycle_biomass_light.png"), dpi=150)
