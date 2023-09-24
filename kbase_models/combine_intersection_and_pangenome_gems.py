import cobra

# Load the models
gem_intersection = cobra.io.read_sbml_model('kbase_models/MIT1002_intersection.sbml')
gem_bayesian_core = cobra.io.load_json_model('kbase_models/alteromonas_bayesian_core.json')

# Label all the metabolites/reactions in the pangenome model as being
# from the pangenome
for met in gem_bayesian_core.metabolites:
    met.annotation['Present in pangenome model'] = True
for rxn in gem_bayesian_core.reactions:
    rxn.annotation['Present in pangenome model'] = True

# Loop through all metabolites in the intersection model
for met in gem_intersection.metabolites:
    # Add to the pangenome model
    if met.id not in gem_bayesian_core.metabolites:
        gem_bayesian_core.add_metabolites(met)
    # Mark as being in the intersection model
    met.annotation['Present in intersection model'] = True

# Loop through all reactions in the intersection model
for rxn in gem_intersection.reactions:
    # Add to the pangenome model
    if rxn.id not in gem_bayesian_core.reactions:
        gem_bayesian_core.add_reaction(rxn)
    # Mark as being in the intersection model
    rxn.annotation['Present in intersection model'] = True

# Give the model a good ID
gem_bayesian_core.id = 'Alteromonas_bayesian_core'

# Save the model file as a JSON file
cobra.io.save_json_model(gem_bayesian_core, 'kbase_models/GEM-MIT1002_intersection_and_pangenome.json')
