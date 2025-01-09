# import os

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


# # Load the biomass compositions from the TSV files
# model_seed_path = "/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase/Templates"

# # Load all of the biomass compositions into one dataframe, and keep only the unique values
# core_template = pd.read_csv(os.path.join(model_seed_path, "Core", "BiomassCompounds.tsv"), sep="\t", header=0)
# gpos_template = pd.read_csv(os.path.join(model_seed_path, "GramPositive", "BiomassCompounds.tsv"), sep="\t", header=0)
# gneg_template = pd.read_csv(os.path.join(model_seed_path, "GramNegative", "BiomassCompounds.tsv"), sep="\t", header=0)
# all_components = pd.concat([core_template, gpos_template, gneg_template])

# Load the producibility results
df = pd.read_csv("biomass/biomass_producibility.csv", header=0, index_col=0)

# Convert the producibility results to a boolean
df_binary = df.replace({True: 1, False: 0})

# Replace NaN with -1
df_binary = df_binary.fillna(-1)

# Define my own color palette
my_palette = sns.color_palette(["gray", "red", "green"])

# Load the model
model = cobra.io.read_sbml_model("2025-01-08_Scott_draft-model-from-KBase.xml")

# Re-name the index to use the names of the metabolites
df_binary.index = [model.metabolites.get_by_id(x + "_c0").name if x + "_c0" in model.metabolites else x for x in df_binary.index]

# Reorder the columns to have gram negative, gram positive, and core in that order
df_binary = df_binary[["gram_negative", "gram_positive", "core"]]

# Plot the producibility results as a heatmap
plt.figure(figsize=(10, 10))
sns.heatmap(df_binary, annot=False, cmap=my_palette, cbar_kws=dict(ticks=[-0.66, 0.0, 0.66], label=["N/A", "Not producible", "Producible"]))

# Add white lines to separate the different cells
for i in range(df_binary.shape[0]):
    plt.axhline(i, color="white", linewidth=0.5)
for i in range(df_binary.shape[1]):
    plt.axvline(i, color="white", linewidth=1)

# Make sure all y ticks/component names are shown
plt.yticks([x + 0.5 for x in range(df_binary.shape[0])], df_binary.index, fontsize=5)

plt.xlabel("Biomass composition")
plt.ylabel("Biomass compound")

plt.tight_layout()

# Save the plot
plt.savefig("biomass/biomass_producibility_heatmap.png")
