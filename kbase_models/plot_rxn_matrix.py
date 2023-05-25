import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load in the reaction matrix
rxn_presence = pd.read_csv('kbase_models/reaction_presence.tsv',
                           header=0,
                           sep='\t')

# Set the reaction column as the index for easy plotting
plotting_matrix = rxn_presence.set_index('Reactions')

# Make a color map
cmap = ['gray', 'green'] # TODO: Pick better colors

# Plot the heatmap
ax = sns.heatmap(plotting_matrix,
                 cmap=cmap,
                 yticklabels=False)

# modify colorbar:
colorbar = ax.collections[0].colorbar
# TODO: Remove magic numbers, these only work when the values are 0 and 1
colorbar.set_ticks([0.25, 0.75])
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
plt.savefig('kbase_models/rxns_in_models.png')