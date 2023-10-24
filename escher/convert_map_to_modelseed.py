import json

# Load in the E. coli Escher map
with open('iJO1366.Central metabolism.json') as f:
    map_json = json.load(f)

# Load in the ModelSEED reactions database
with open('../../ModelSEEDDatabase/Biochemistry/reactions.json') as f:
    reactions = json.load(f)

# Load in the ModelSEED compounds database
with open('../../ModelSEEDDatabase/Biochemistry/compounds.json') as f:
    compounds = json.load(f)

# Loop through the reactions in the map and add the ModelSEED IDs
for reaction in map_json[1]['reactions']:
    # Get the reaction ID
    bigg_id = reaction['bigg_id']

    # Get the ModelSEED ID
    modelseed_id = reaction['bigg_id']

    # Get the ModelSEED reaction
    modelseed_reaction = reactions[modelseed_id]

    # Get the ModelSEED compounds
    modelseed_compounds = modelseed_reaction['compound_ids']

    # Loop through the compounds in the reaction
    for compound in reaction['metabolites']:
        # Get the compound ID
        compound_id = compound['bigg_id']

        # Get the ModelSEED ID
        modelseed_compound_id = modelseed_compounds[compound_id]

        # Get the ModelSEED compound
        modelseed_compound = compounds[modelseed_compound_id]

        # Add the ModelSEED ID to the compound
        compound['bigg_id'] = modelseed_compound_id

        # Add the ModelSEED name to the compound
        compound['name'] = modelseed_compound['name']

        # Add the ModelSEED formula to the compound
        compound['formula'] = modelseed_compound['formula']

        # Add the ModelSEED charge to the compound
        compound['charge'] = modelseed_compound['charge']