# Make an escher map file with all of the reactions in the model using random
# coordinates
# Also add text fields for human-readable reaction and metabolite names

import json

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

# For a reaction
# Are there multiple reactants
    # If there is only one reactant, can go straight from metabolite to midmarker
    # If there are multiple reactants, need to go from metabolite to multimarker to midmarker
# Are there multiple products
    # If there is only one product, can go straight from midmarker to metabolite
    # If there are multiple products, need to go from midmarker to multimarker to metabolite

# Save the escher map json file
with open('escher/random_map.json', 'w') as f:
    json.dump(map_file, f, indent=4)