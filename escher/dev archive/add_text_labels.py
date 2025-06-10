import json

# Read in the map file, assume that the human readable names are already
# in the "name" field
with open('escher/e_coli.Core metabolism.modelseed.json') as f:
    map_file = json.load(f)

# Set if you want to add text labels to the reactions (as well as the metabolites)
rxn_labels = False

# Function
def add_text_label(map_file, id, x, y, name):
    # Add a text label with the same name and the same x and y coordinates
    map_file[1]['text_labels'][id] = {
        "text": name,
        "x": x,
        "y": y
    }

# Get the current max ID in the whole map file
max_id = int(max(list(map_file[1]['nodes'].keys()) +
              list(map_file[1]['reactions'].keys()) +
              list(map_file[1]['text_labels'].keys())))

# Loop through all nodes
for node in map_file[1]['nodes']:
    # Get the node info
    node_info = map_file[1]['nodes'][node]
    # If node doesn't have a "name" field, skip
    if 'name' not in node_info:
        continue
    # If node does have a "name" field, add a text label with the same name
    # and the same x and y coordinates
    add_text_label(map_file,
                   str(max_id + 1),
                   node_info['label_x'],
                   node_info['label_y'] - 20, # So the text label is above the ID
                   node_info['name'])
    max_id += 1

# Loop through all the reactions (if rxn_labels is True)
if rxn_labels:
    for reaction in map_file[1]['reactions']:
        rxn_info = map_file[1]['reactions'][reaction]
        # All reactions should have a name, but maybe good to check
        # Add a text label with the same name and the same x and y coordinates
        add_text_label(map_file,
                    str(max_id + 1),
                    rxn_info['label_x'],
                    rxn_info['label_y'] - 20,
                    rxn_info['name'])
        max_id += 1

# Write out the map file
with open('escher/e_coli_core.Coremetabolism.modelseed.With labels.json', 'w') as f:
    json.dump(map_file, f, indent=4)
