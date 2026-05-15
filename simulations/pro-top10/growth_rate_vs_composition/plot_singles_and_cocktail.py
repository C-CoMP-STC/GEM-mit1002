from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set file paths
FILE_PATH = Path(__file__).resolve().parent
OUT_PATH = FILE_PATH / "results"
TOP_10_DIR = FILE_PATH.parent
REPO_ROOT = FILE_PATH.parents[2]
TEST_FILE_DIR = REPO_ROOT / "test" / "test_files"

# Make the results directory if it doesn't exist
OUT_PATH.mkdir(exist_ok=True)

# Load the results from the single substrate and cocktail simulations
results = pd.read_csv(OUT_PATH / "single_substrate_results.csv")

# Plot the results
plt.figure(figsize=(8, 6))
sns.barplot(x="substrate", y="growth_rate", data=results)
plt.xlabel("Substrate")
plt.ylabel("Growth Rate (1/hr)")
plt.title("Growth Rate on Single Substrates and Cocktail")
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(OUT_PATH / "growth_rate_comparison.png")
