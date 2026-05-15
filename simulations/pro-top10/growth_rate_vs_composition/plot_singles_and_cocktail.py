from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent
REPO_ROOT = FILE_PATH.parents[2]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"

# Add the repo root to the system path so we can import from the plot_styles file
sys.path.append(str(REPO_ROOT))
import plot_styles  # Import the plot styles from the repo

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)

# Load the results from the single substrate and cocktail simulations
results = pd.read_csv(OUT_PATH / "single_substrate_results_total_c_60.csv")

# Select the colors I want to use from the color palette
colors = [
    plot_styles.tsitp_colors["dark_blue"],
    plot_styles.tsitp_colors["dark_green"],
    plot_styles.tsitp_colors["dark_orange"],
]

# Plot the results
fig, ax = plt.subplots(figsize=(8, 6))
sns.barplot(
    x="substrate",
    y="growth_rate",
    data=results,
    hue="condition",
    ax=ax,
    palette=colors,
)
plt.xlabel("Substrate")
plt.ylabel("Growth Rate (1/hr)")
plt.title("Growth Rate on Single Substrates and Cocktail")
plt.xticks(rotation=45)

# Move the legend to the right of the plot
ax.legend(loc="upper center", bbox_to_anchor=(1.15, 0.5), ncol=1, frameon=False)

plt.tight_layout()

# Apply the plot style
plot_styles.set_plot_style(ax)

plt.savefig(OUT_PATH / "growth_rate_comparison.png")
