import cobra
import pandas as pd 

# Load the models
gem_dram = cobra.io.read_sbml_model('kbase_models/GEM-MIT1002_DRAM.sbml')
gem_prokka = cobra.io.read_sbml_model('kbase_models/GEM-MIT1002_Prokka.sbml')
gem_rast = cobra.io.read_sbml_model('kbase_models/GEM-MIT1002_RAST.sbml')

gem_union = cobra.io.read_sbml_model('kbase_models/MIT1002_union.sbml')
gem_intersection = cobra.io.read_sbml_model('kbase_models/MIT1002_intersection.sbml')

gem_2_or_more = cobra.io.read_sbml_model('kbase_models/MIT1002_2OrMore.sbml')
gem_majority_rule = cobra.io.read_sbml_model('kbase_models/MIT1002_MajorityRule.sbml')
gem_priority_list = cobra.io.read_sbml_model('kbase_models/MIT1002_PriorityList.sbml')

gem_bayesian_core = cobra.io.load_json_model('kbase_models/alteromonas_bayesian_core.json')

models_with_names = {'DRAM': gem_dram,
                     'Prokka': gem_prokka,
                     'RASTtk': gem_dram,
                     'Intersection': gem_intersection,
                     '2+ Annotations': gem_2_or_more,
                     'Priority List': gem_priority_list,
                     'Majority Rule': gem_majority_rule,
                     'Union': gem_union,
                     'Bayesian Pangenome Core': gem_bayesian_core
                     }

# Load the ModelSEED Pathways table
pathways_db = pd.read_csv('../../ModelSEEDDatabase/Biochemistry/Pathways/ModelSEED_Subsystems.tsv',
                          header=0,
                          sep='\t')

# Filter to only the reactions in the union model
# FIXME: To get the list of all possible reactions, combine the list from the
# union model and the pangenome
all_rxns = list(set([r.id[0:-3] for r in gem_union.reactions]) | set([r.id for r in gem_bayesian_core.reactions]))
pathways_db = pathways_db[pathways_db['Reaction'].isin(all_rxns)]

# Add if the reaction is present in the different models to the table
for name in models_with_names:
    model = models_with_names[name]
    rxn_presence = []
    # Quick and dirty way to deal with the different nomenclature, just handle
    # the Bayesian core model totally separated
    if name == 'Bayesian Pangenome Core':
        for index, row in pathways_db.iterrows():
            if row['Reaction'] in [r.id for r in model.reactions]:
                rxn_presence.append(1)
            else:
                rxn_presence.append(0)
    else:
        for index, row in pathways_db.iterrows():
            if row['Reaction'] + '_c0' in [r.id for r in model.reactions]: # FIXME: The _c0 is what is making it think that there are no reactions in the bayesian core
                rxn_presence.append(1)
            else:
                rxn_presence.append(0)
    pathways_db[name] = rxn_presence

# Save as a CSV
pathways_db.to_csv('kbase_models/pathway_reaction_presence.tsv',
                   sep='\t',
                   index=False)