import os
import sys

import cobra
import pandas as pd

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

# Load the results
results = pd.read_csv(os.path.join(OUT_DIR, "results.csv"))

# Stacked bar plot of the carbon fates for the different conditions
data = results.set_index("sim_name")[["biomass", "organic_c", "co2"]]
g = carbon_fates_bar(data)
plt.savefig(os.path.join(output_folder, "carbon_fates.png"))

# Stacked bar plot of the normalized carbon fates for the different conditions
data_norm = results.set_index("sim_name")[["biomass_norm", "organic_c_norm", "co2_norm"]]
# Rename the columns to match the function
data_norm.columns = ["biomass", "organic_c", "co2"]
g_norm = carbon_fates_bar(data_norm)
plt.savefig(os.path.join(output_folder, "carbon_fates_norm.png"))

# Subset the results to only include the pFBA on realistic O2
clean_data = data_norm[data_norm.index.str.contains("real_o2_pfba")]
g_clean = carbon_fates_bar(clean_data)
# Relabel the x tick labels
g_clean.set_xticklabels(
    ["Glucose", "Acetate", "Heavy Glucose Mix", "Heavy Acetate Mix"]
)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "carbon_fates_norm_clean.png"))

