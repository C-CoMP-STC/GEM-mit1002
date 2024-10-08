import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the known growth phenotypes with predicted growth phenotypes
growth_phenotypes = pd.read_csv('kbase_models/known_growth_phenotypes_w_pred.tsv',
                                sep='\t')

# Clean up the data frame
growth_phenotypes = growth_phenotypes.set_index('c_source')
growth_phenotypes = growth_phenotypes.drop(columns=['met_id', 'notes'])

# Rename the columns and the index to be longer/more descriptive
growth_phenotypes.index.name = 'Carbon Source'
growth_phenotypes = growth_phenotypes.rename(columns={'exp_growth': 'Experimental Growth'})

# Replace all of the "No Exchange" values with "No"
growth_phenotypes = growth_phenotypes.replace('No Exchange', 'No')

# Make a dictionary for the phenotypes to numbers
value_to_int = {'Unsure': 0, 'No': 1, 'Yes': 2}
n = len(value_to_int)

# Make a colormap of specified colors (in numerical order for the phenotypes)
# cmap = ['gray', '#F18F01', '#399E5A'] # Gray, orange, green
cmap = ['#5E5E5E', '#FF7D0A', '#024064'] # C-CoMP gray, orange, and dark blue

# Plot the heatmap
ax = sns.heatmap(growth_phenotypes.replace(value_to_int),
                 cmap=cmap,
                 linewidths=4,
                 linecolor='white')

# modify colorbar:
colorbar = ax.collections[0].colorbar 
r = colorbar.vmax - colorbar.vmin 
colorbar.set_ticks([colorbar.vmin + r / n * (0.5 + i) for i in range(n)])
colorbar.set_ticklabels(list(value_to_int.keys()))

# Move the x-axis labels to the top
plt.tick_params(axis='both',
                which='major',
                labelsize=10,
                labelbottom = False,
                bottom=False,
                top = False,
                labeltop=True)
plt.tick_params(axis='x',
                labelrotation = 90)

# Make sure that the y-axis labels are not cut off
plt.tight_layout()

# Save the figure
plt.savefig('kbase_models/exp_vs_pred_growth_phenotypes.png')
