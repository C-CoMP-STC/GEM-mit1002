import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the known growth phenotypes with predicted growth phenotypes
growth_phenotypes = pd.read_csv('known_growth_phenotypes_w_pred.tsv', sep='\t')

# Plot a categorical heatmap of the growth phenotypes, where the rows
# are the metabolites and the columns are the experimental and predicted
# growth phenotypes. Show growth as green and no growth as red, unsure
# as yellow, and no exchange reaction as gray.
# First, make a new dataframe with the metabolites as the rows and the
# experimental and predicted growth phenotypes as the columns
growth_phenotypes = growth_phenotypes.set_index('c_source')
growth_phenotypes = growth_phenotypes[['exp_growth', 'pred_growth']]

# Make a dictionary for the phenotypes to numbers
value_to_int = {'No': 0, 'Yes': 1, 'Unsure': 2, 'No Exchange': 3}
n = len(value_to_int)

# Make a colormap of specified colors
cmap = ['red', 'green', 'yellow', 'gray']

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

# Save the figure
plt.savefig('exp_vs_pred_growth_phenotypes.png')
