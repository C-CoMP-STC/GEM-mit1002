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

# Plot results
# Identify media metabolites with meaningful concentration changes
media["time"] = media["cycle"] * TIME_STEP
met_range = media.groupby("metabolite")["conc_mmol"].apply(lambda x: x.max() - x.min())
active_mets = met_range[met_range > 1e-10].index.tolist()

n_panels = 1 + len(active_mets)
fig, axes = plt.subplots(n_panels, 1, figsize=(10, 3 * n_panels), sharex=True)
if n_panels == 1:
    axes = [axes]

time_hours = biomass["cycle"] * TIME_STEP
axes[0].plot(time_hours, biomass.iloc[:, 1])
axes[0].set_ylabel("Biomass (gDW)")
axes[0].set_title("Prochlorococcus diel cycle simulation")

for ax, met in zip(axes[1:], active_mets):
    met_data = media[media["metabolite"] == met]
    ax.plot(met_data["time"], met_data["conc_mmol"])
    ax.set_ylabel("mmol")
    ax.set_title(f"{met_names.get(met, met)} ({met})")

axes[-1].set_xlabel("Time (hours)")
plt.tight_layout()
plt.savefig(os.path.join(FILE_DIR, "diel_cycle_results.png"), dpi=150)
