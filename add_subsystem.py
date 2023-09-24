import json

# Open the model json file
with open('model.json', 'r') as f:
    model = json.load(f)

# Loop through all of the reactions
for reaction in model['reactions']:
    reaction['subsystem'] = "TBD"

# Save the model file
with open('model.json', 'w') as f:
    json.dump(model, f, indent=4)

# Once the subsystem field has been added, fill it in
