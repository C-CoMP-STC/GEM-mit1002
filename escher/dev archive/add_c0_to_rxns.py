import json

# Load the map
with open('escher/e_coli.Core metabolism.modelseed.json') as f:
    map = json.load(f)

# Loop through all the reactions
for reaction in map[1]['reactions']:
    rxn_info = map[1]['reactions'][reaction]
    # Add "_c0" to the reaction ID
    map[1]['reactions'][reaction]['bigg_id'] += '_c0'

# Save the new map
with open('escher/e_coli.Core metabolism.modelseed.json', 'w') as f:
    json.dump(map, f, indent=4)
