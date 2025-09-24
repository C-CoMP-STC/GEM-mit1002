import cobra
import json

# Load the model
model = cobra.io.read_sbml_model('model.xml')

# Load the map json file
with open('iJO1366.Central metabolism.json') as f:
    map_json = json.load(f)

# Make a new empty json object for the ALT map, with all of the necessary
# header information
# Not sure what to put for homepage, so I ommitted it
alt_map_json = [{'map_name': 'mit1002.Central metabolism',
                 'map_id': '',
                 'map_description': 'Central metabolism (MIT1002)',
                 'schema': 'https://escher.github.io/escher/jsonschema/1-0-0#'
                },
                {'reactions': {},
                 'nodes': {},
                 'text_labels': {},
                 'canvas': {}
                }]

# Loop through all the reactions in the map
for rxn in map_json[1]['reactions']:
    # Get the BiGG ID of the reaction
    bigg_id = map_json[1]['reactions'][rxn]["bigg_id"]
    # If the reaction is in the model, add it to a new json file
    try:
        model_rxn = model.reactions.get_by_id(bigg_id)
        # print(bigg_id + " is in the model")
        # Add it to the new json file
        alt_map_json[1]['reactions'][rxn] = map_json[1]['reactions'][rxn]
    except:
        # print(bigg_id + " is NOT in the model")
        continue

# Loop through all the nodes in the map
for node in map_json[1]['nodes']:
    # If the node is a metabolite node, check that it is in the model
    if map_json[1]['nodes'][node]["node_type"] == "metabolite":
        # Get the BiGG ID of the metabolite
        bigg_id = map_json[1]['nodes'][node]["bigg_id"]
        # If the reaction is in the model, add it to a new json file
        try:
            model_rxn = model.metabolites.get_by_id(bigg_id)
            # print(bigg_id + " is in the model")
            # Add it to the new json file
            alt_map_json[1]['nodes'][node] = map_json[1]['nodes'][node]
        except:
            # print(bigg_id + " is NOT in the model")
            continue
    # Otherwise (i.e. for multimarkers or midmarkers), add it to the new json file
    else:
        alt_map_json[1]['nodes'][node] = map_json[1]['nodes'][node]

# Copy over the text labels and canvas exactly as thet are in the original map
# TODO: Remove the unnecessary text labels
alt_map_json[1]['text_labels'] = map_json[1]['text_labels']
alt_map_json[1]['canvas'] = map_json[1]['canvas']

# Serializing json
json_object = json.dumps(alt_map_json, indent=4)
 
# Writing to sample.json
with open("alt_map.json", "w") as outfile:
    outfile.write(json_object)