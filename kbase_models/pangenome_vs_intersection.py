import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import cobra

# Load the models
gem_intersection = cobra.io.read_sbml_model('kbase_models/MIT1002_intersection.sbml')
gem_bayesian_core = cobra.io.load_json_model('kbase_models/alteromonas_bayesian_core.json')

# Make a file to list the reactions that the two models have in common
rxn_file = open("kbase_models/pangenome_vs_intersection_rxns.txt", "w")

# Get the reactions that the two have in common
for reaction in gem_intersection.reactions:
    if reaction.id[0:-3] in [r.id for r in gem_bayesian_core.reactions]:
        if 'bigg.reaction' in reaction.annotation.keys():
            bigg_id = reaction.annotation['bigg.reaction']
        else:
            bigg_id = 'No BiGG ID'
        lower_bound = reaction.lower_bound
        upper_bound = reaction.upper_bound
        readable_reaction_string = " + ".join([f"{coefficient} {met.name}" for met, coefficient in reaction.metabolites.items()])
        rxn_file.write(f"{reaction.id} | {bigg_id} | Bounds: [{lower_bound}, {upper_bound}] | {reaction.name} | {readable_reaction_string}\n")

# Close the file
rxn_file.close()

# Make a venn diagram
venn2([set([r.id[0:-3] for r in gem_intersection.reactions]),
       set([r.id for r in gem_bayesian_core.reactions])],
       set_labels=('Intersection of All \n KBase Models',
                   'Pangenome \n Bayesian \n Core Model')
       )

# Style
plt.tight_layout()  # This alone doesn't keep labels from hanging off

# Save the figure
plt.savefig('kbase_models/pangenome_vs_intersection.png')