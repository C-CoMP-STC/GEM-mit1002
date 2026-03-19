import os
import sys

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))

# Add the project root to the system path to import my plot_styles file
sys.path.append(PROJECT_ROOT)
# Import the plot style from the plot_styles.py file
from plot_styles import ccomp_colors

# Load the data
data = pd.read_csv(os.path.join(FILE_DIR, "growth_match_summary.csv"))

# Skip rows that say "Error" in the % Match column
data = data[data["% Match"] != "ERROR"]

# Convert columns to numeric type
# TODO: Do I need any other columns numerica? Are there any that can't be?
cols_to_fix = ["Matches", "% Match", "PR Number", "Unbounded Flux Reactions"]
data[cols_to_fix] = data[cols_to_fix].apply(pd.to_numeric, errors='coerce')

# Sort the data by PR number and reset the index
data = data.sort_values("PR Number").reset_index(drop=True)

# Create a figure with twin y axes
fig, ax1 = plt.subplots(figsize=(10, 6))
ax2 = ax1.twinx()
# Scale the unbounded flux axis to show detail near zero
ax2.set_yscale('symlog', linthresh=1)

# Plot the number of matches on the left axis
ax1.plot(
    data.index,
    data["Matches"],
    marker="o",
    linestyle="-",
    color=ccomp_colors["dark_blue"],
)

# Plot the number of arbitrarily large reactions on the right axis
ax2.plot(
    data.index,
    data["Unbounded Flux Reactions"],
    marker="o",
    linestyle="-",
    color=ccomp_colors["orange"],
)

# Make it so it looks like the matches line is "on top"
# Move ax1 to a higher z-order than ax2
ax1.set_zorder(ax2.get_zorder() + 1)
# Make ax1's background transparent so ax2 is still visible behind it
ax1.patch.set_visible(False)

# Titles
plt.title("Model Performace Over Time")
plt.xlabel("Pull Request Number")
# Show x-ticks for every point, but label them with the actual PR number not the index in the df
plt.xticks(data.index, data["PR Number"], rotation=90)
# Set y axes labels
ax1.set_ylabel("Number of Growth Phenotypes Matching Experimental Data")
ax1.set_ylim(0, 50)  # See the full range
ax2.set_ylabel("Number of Unique Reactions with Flux > 100 (Log Scale)")

# Save the plot
plt.savefig(
    os.path.join(FILE_DIR, "match_over_time.png"),
)
