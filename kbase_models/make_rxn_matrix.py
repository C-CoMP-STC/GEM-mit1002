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

# Load the model seed database for the gram negative bacteria template
gram_neg_db = pd.read_csv('../../ModelSEEDDatabase/Templates/GramNegative/Reactions.tsv',
                          sep='\t',
                          header=0)

# Filter the bayesian core model to only include reactions that are in the
# gram negative database
rxns_to_remove = [r.id for r in gem_bayesian_core.reactions if
                  r.id not in gram_neg_db['id'].tolist() and
                  not r.id.startswith('EX_')] # Don't think this was necessary, since it didn't have any exchange reactions to begin with
gem_bayesian_core.remove_reactions(rxns_to_remove)

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

# Get a list of all the possible reactions in all of the models
# ( Really just the combon of the union model and the pangenome)
all_rxns = list(set([r.id[0:-3] for r in gem_union.reactions if not r.id.endswith("_b")]) | set([r.id for r in gem_bayesian_core.reactions]))

df = pd.DataFrame({'Reactions': all_rxns})

# Add if the reaction is present in the different models to the table
for name in models_with_names:
    model = models_with_names[name]
    rxn_presence = []
    # Quick and dirty way to deal with the different nomenclature, just handle
    # the Bayesian core model totally separated
    if name == 'Bayesian Pangenome Core':
        for rxn in all_rxns:
            if rxn in [r.id for r in model.reactions]:
                rxn_presence.append(1)
            else:
                rxn_presence.append(0)
    else:
        for rxn in all_rxns:
            if rxn + '_c0' in [r.id for r in model.reactions]:
                if rxn.startswith('EX_'):
                    rxn_presence.append(2)
                else:
                    rxn_presence.append(1)
            elif rxn + '_e0' in [r.id for r in model.reactions]:
                if rxn.startswith('EX_'):
                    rxn_presence.append(2)
                else:
                    rxn_presence.append(1)
            elif rxn in [r.id for r in model.reactions]:
                if rxn.startswith('EX_'):
                    rxn_presence.append(2)
                else:
                    rxn_presence.append(1)
            else:
                rxn_presence.append(0)
    df[name] = rxn_presence

# Save as a CSV
df.to_csv('kbase_models/reaction_presence.tsv',
          sep='\t',
          index=False)