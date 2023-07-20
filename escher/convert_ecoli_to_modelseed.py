import json

# Load the iJO1366.Central metabolism map
with open('escher/iJO1366.Central metabolism.json') as f:
    map = json.load(f)

# Load the ModelSeed compound database
with open('../../ModelSEEDDatabase/Biochemistry/compounds.json') as f:
    compounds = json.load(f)

# Load the ModelSeed reaction database
with open('../../ModelSEEDDatabase/Biochemistry/reactions.json') as f:
    reactions = json.load(f)

# Helper functions
# Convert the alias string into a dictionary
def parse_aliases(alias_list):
    alias_dict = {}
    for meta_list in alias_list:
        # Stupid way to get around incase there are colons in the alias list
        meta_name = meta_list.split(':')[0]
        meta_entries = ':'.join(meta_list.split(':')[1:])
        # Seperate the alias IDs by semicolons and remove whitespace
        alias_ids = [x.strip() for x in meta_entries.split(';')]
        # Add the alias IDs to the dictionary
        alias_dict[meta_name] = alias_ids
    return alias_dict

# For every reaction in the map
for reaction in map[1]['reactions']:
    rxn_info = map[1]['reactions'][reaction]
    bigg_id = rxn_info['bigg_id']
    # Search the ModelSeed reaction database for the reaction
    potential_rxns = []
    for rxn in reactions:
        # Convert the alias string into a dictionary
        if 'aliases' not in rxn.keys() or rxn['aliases'] is None:
            continue
        alias_dict = parse_aliases(rxn['aliases'])
        # If the bigg_id is in the alias dictionary, add it to the list
        if 'BiGG' in alias_dict.keys() and bigg_id in alias_dict['BiGG']:
            potential_rxns.append(rxn)
    # If there is only one match, use it
    if len(potential_rxns) == 1:
        print("Success: " + bigg_id)
        modelseed_rxn = potential_rxns[0]
    # If there are multiple matches, give a warning and skip it
    elif len(potential_rxns) > 1:
        print("WARNING: Multiple matches for " + bigg_id)
        continue
    # Change the reaction name (bigg_id) to the ModelSeed ID
    rxn_info['bigg_id'] = modelseed_rxn['id']

# For every node in the map
for node in map[1]['nodes']:
    node_info = map[1]['nodes'][node]
    # Use the "node type" field to determine if it is a metabolite or a
    # reaction
    if node_info['node_type'] == 'metabolite':
        # Remove the compartment from the metabolite name
        # Hardcoding assuming that all compartments are 1 character long
        # and are preceded by an underscore (i.e. _c, _e)
        bigg_id = node_info['bigg_id'][:-2]
        compartment = node_info['bigg_id'][-1]
        # Search the ModelSeed compound database for the metabolite
        potential_metabs = []
        for met in compounds:
            # Convert the alias string into a dictionary
            if 'aliases' not in met.keys() or met['aliases'] is None:
                continue
            alias_dict = parse_aliases(met['aliases'])
            # If the bigg_id is in the alias dictionary, add it to the list
            if 'BiGG' in alias_dict.keys() and bigg_id in alias_dict['BiGG']:
                potential_metabs.append(met)
        # If there is only one match, use it
        if len(potential_metabs) == 1:
            modelseed_met = potential_metabs[0]
        # If there are multiple matches, give a warning and skip it
        elif len(potential_metabs) > 1:
            print("WARNING: Multiple matches for " + bigg_id)
            continue
        # Add the compartement information back to the ID, but in the new style
        # (i.e. _c --> _c0)
        if compartment == 'c':
            modelseed_id = modelseed_met['id'] + '_c0'
        elif compartment == 'e':
            modelseed_id = modelseed_met['id'] + '_e0'
        else:
            print("WARNING: Unknown compartment " + compartment)
            continue
    else:
        # If it's something else (i.e. a multimarker), skip it
        continue
    # Change the node name (bigg_id) to the ModelSeed ID
    # FIXME: Can I call the ModelSeed ID something else so it isn't confusing?
    node_info['bigg_id'] = modelseed_id

# Save the new map
with open('escher/iJO1366.Central metabolism.modelseed.json', 'w') as f:
    json.dump(map, f, indent=2)