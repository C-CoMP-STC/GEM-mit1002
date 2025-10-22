import os

import matplotlib.pyplot as plt
import pandas as pd

FILE_DIR = os.path.dirname(os.path.abspath(__file__))

# Load the data
data = pd.read_csv(os.path.join(FILE_DIR, "growth_match_summary.csv"))

# Removing this for when showing the Acetate/Leucine/Isoleucine problem since
# a lot of Devlin's recent changes are for GPRs
# # Loop through the rows, if the row does not have a different number of
# # reactions, metabolites, or genes (i.e. if the model has not changed),
# # remove it
# data = data[
#     (data["Reactions"] != data["Reactions"].shift(1))
#     | (data["Metabolites"] != data["Metabolites"].shift(1))
#     | (data["Genes"] != data["Genes"].shift(1))
# ]

# Skip rows that say "Error" in the % Match column
data = data[data["% Match"] != "ERROR"]

# Convert the % Match and the PR number columns to numeric
data["% Match"] = pd.to_numeric(data["% Match"])
data["PR Number"] = pd.to_numeric(data["PR Number"])

# Plot the percentage of matches
plt.figure(figsize=(10, 6))
plt.plot(
    data["PR Number"],
    data["% Match"],
    marker="o",
    linestyle="-",
    color="blue",
)
plt.title("Percentage of Matches Over Time")
plt.xlabel("Pull Request Number")
# Make sure x ticks are integers
plt.xticks(data["PR Number"].unique())
plt.ylabel("Percentage Match (%)")

# Save the plot
plt.savefig(
    os.path.join(FILE_DIR, "match_over_time.png"),
)
