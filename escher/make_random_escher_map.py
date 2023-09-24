# Make an escher map file with all of the reactions in the model using random
# coordinates
# Also add text fields for human-readable reaction and metabolite names

import json
import cobra

# Make the start of the escher map dictionary
map_file = []
map_file.append({
        "homepage": "https://escher.github.io",
        "map_description": "Randomly generated locations for all reactions in the ALT model",
        "map_id": "",
        "map_name": "MIT1002.Random map",
        "schema": "https://escher.github.io/escher/jsonschema/1-0-0#"
    })
map_file.append({"canvas": {
            "height": 6642.899758,
            "width": 9047.498999,
            "x": -286.3920445,
            "y": 132.78023100000001
        },
        "nodes": {},
        "reactions": {},
        "text_labels": {}})

# Load the ALT model
with open('model.json') as f:
    model = json.load(f)

# Start a node counter
node_counter = 1

# For a reaction
for reaction in model['reactions']:
    # Determine the number of reactants and products using the sign of the
    # stoichiometric coefficients (negative is a reactant, positive is a
    # product)
    reactants = [met for met in reaction['metabolites'] if reaction['metabolites'][met] < 0]
    products = [met for met in reaction['metabolites'] if reaction['metabolites'][met] > 0]
    # Determine the number of interior segments needed for the reaction
    # arrow.
    n_segments = 1 # Every reaction has to have a mid marker
    if len(reactants) > 1:
        # If there are multiple reactants, need a multimarker to connect
        # all of the reactants to the mid-marker
        n_segments += 1 
    if len(products) > 1:
        # If there are multiple products, need a multimarker to connect
        # all of the products to the mid-marker
        n_segments += 1
    # TODO: Stop hardocoding the starting positions
    starting_x = 500
    starting_y = 500
    # Determine the center of the reaction
    reactants_x = starting_x + len(reactants) * 10
    products_x = starting_x + len(products) * 10
    if len(reactants) > len(products):
        
    # Make a node for each reactant
    for idx, met in enumerate(reactants):
        # Get the metabolite info
        met_info = [met for met in model['metabolites'] if met['id'] == met]
        map_file[1]['nodes'][met] = { # Not sure if the node id should be the metabolite id or the node counter
                "bigg_id": met,
                "label_x": starting_x + idx * 10,
                "label_y": starting_y + 10,
                "name": model.metabolites[met]['name'],
                "node_is_primary": True,
                "node_type": "metabolite",
                "x": starting_x + idx * 10,
                "y": starting_y
            }
    # Get the center of the all the reactant nodes

    # Determine the position for the nodes
    # Make the reaction dictionary


# Save the escher map json file
with open('escher/random_map.json', 'w') as f:
    json.dump(map_file, f, indent=4)