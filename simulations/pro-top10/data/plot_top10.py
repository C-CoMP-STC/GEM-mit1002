import os

import matplotlib.pyplot as plt
import pandas as pd

FILE_DIR = os.path.dirname(__file__)
REPO_DIR = os.path.dirname(os.path.dirname(FILE_DIR))
OUT_DIR = FILE_DIR
# Make the output directory if it doesn't exist
os.makedirs(OUT_DIR, exist_ok=True)

# Load the top 10 metabolite file
top_10_metabolites = pd.read_csv(os.path.join(FILE_DIR, "top10_exometabolites.csv"))

# Filter the top 10 metabolites to only include those from Prochlorococcus marinus
top_10_metabolites = top_10_metabolites[
    top_10_metabolites["organism"] == "Prochlorococcus marinus"
].copy()

# Sort the dataframe by the metabolite concentration
# Descending so that in the bar chart the largest metabolite will be on the left
top_10_metabolites.sort_values("carbon_concentration", ascending=False, inplace=True)

# Transpose the data so that I have the metabolites as columns
df_plot = top_10_metabolites.set_index("metabolite")[["carbon_concentration"]].T

# Plot a horizontal stacked bar graph of the top 10 metabolites of their metabolite concentations
# Just one bar for the whole organism, the width of the segment will correspond to the concentration of that metabolite, and the color will correspond to the metabolite class
df_plot.plot(kind="barh", stacked=True, legend=True, figsize=(10, 2))

# move the legend outside of the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

# Save
plt.savefig(os.path.join(OUT_DIR, "top10_metabolites.png"), bbox_inches="tight")
