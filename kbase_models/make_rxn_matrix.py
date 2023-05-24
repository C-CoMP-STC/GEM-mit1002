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

models_with_names = {'DRAM': gem_dram,
                     'Prokka': gem_prokka,
                     'RASTtk': gem_dram,
                     'Intersection': gem_intersection,
                     '2+ Annotations': gem_2_or_more,
                     'Priority List': gem_priority_list,
                     'Majority Rule': gem_majority_rule,
                     'Union': gem_union,
                     }

# Load the ModelSEED Pathways table
pathways_db = pd.read_csv('../../ModelSEEDDatabase/Biochemistry/Pathways/ModelSEED_Subsystems.tsv',
                          header=0,
                          sep='\t')

# Filter to only the reactions in the union model
pathways_db = pathways_db[pathways_db['Reaction'].isin([r.id[0:-3] for r in gem_union.reactions])]

# Add if the reaction is present in the different models to the table
for name in models_with_names:
    model = models_with_names[name]
    rxn_presence = []
    for index, row in pathways_db.iterrows():
        if row['Reaction'] + '_c0' in [r.id for r in model.reactions]:
            rxn_presence.append(1)
        else:
            rxn_presence.append(0)
    pathways_db[name] = rxn_presence

# Save as a CSV
pathways_db.to_csv('kbase_models/pathway_reaction_presence.tsv',
                   sep='\t',
                   index=False)