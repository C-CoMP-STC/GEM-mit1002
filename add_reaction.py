import cobra
from gem_utilities import curation
import json


# Load the current version of the model
old_model = cobra.io.read_sbml_model('model.xml')

# Load the ModelSEED databases
rxn_db = json.load(
    open('/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase/Biochemistry/reactions.json')
)
met_db = json.load(
    open('/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase/Biochemistry/compounds.json')
)

# Convert the ModelSEED databases to dictionaries for easy searching
template_rxn_db = {rxn["id"]: rxn for rxn in rxn_db if not rxn["is_obsolete"]}
template_met_db = {met["id"]: met for met in met_db if not met["is_obsolete"]}

new_model = curation.add_ms_reaction_from_id(old_model, template_rxn_db, template_met_db, 'rxn00333')

# Save the new model
cobra.io.write_sbml_model(new_model, 'model.xml')
