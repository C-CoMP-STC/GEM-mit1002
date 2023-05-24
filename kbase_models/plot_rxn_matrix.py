import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load in the reaction matrix
rxn_presence = pd.read_csv('kbase_models/pathway_reaction_presence.tsv',
                           header=0,
                           sep='\t')

# Set the reaction matrix as the index and subset just the model columns
# For easy plotting
plotting_matrix = rxn_presence.set_index('Reaction')
plotting_matrix = plotting_matrix[list(rxn_presence.columns)[5:]]

# Make a color map
cmap = ['gray', 'green'] # TODO: Pick better colors

# Plot the heatmap
ax = sns.heatmap(plotting_matrix,
                 cmap=cmap)

# modify colorbar:
colorbar = ax.collections[0].colorbar 
colorbar.set_ticks([0.25, 0.75]) # TODO: Remove magic numbers, these only work when the values are 0 and 1
colorbar.set_ticklabels(['No', 'Yes'])

# Move the x-axis labels to the top
plt.tick_params(axis='both',
                which='major',
                labelsize=10,
                labelbottom = False,
                bottom=False,
                top = False,
                labeltop=True,
                labelrotation = 90)

# Make sure that the y-axis labels are not cut off
plt.tight_layout()

# Save the figure
plt.savefig('kbase_models/rxns_in_models_by_pathway.png')