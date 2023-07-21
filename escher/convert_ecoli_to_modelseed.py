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

# Make a dictionary for all the compounds I change the ID of
bigg_to_modelseed = {}

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
            modelseed_id = modelseed_met['id'] + '_' + compartment
    else:
        # If it's something else (i.e. a multimarker), skip it
        continue
    # Add the old and new IDs to the dictionary
    bigg_to_modelseed[node_info['bigg_id']] = modelseed_id
    # Change the node name (bigg_id) to the ModelSeed ID
    # FIXME: Can I call the ModelSeed ID something else so it isn't confusing?
    node_info['bigg_id'] = modelseed_id

# For every reaction in the map
for reaction in map[1]['reactions']:
    rxn_info = map[1]['reactions'][reaction]
    bigg_id = rxn_info['bigg_id']
    # Search the ModelSeed reaction database for the reaction
    potential_rxns = {}
    for rxn in reactions:
        # Convert the alias string into a dictionary
        if 'aliases' not in rxn.keys() or rxn['aliases'] is None:
            continue
        alias_dict = parse_aliases(rxn['aliases'])
        # If the bigg_id is in the alias dictionary, add it to the list
        if 'BiGG' in alias_dict.keys() and bigg_id in alias_dict['BiGG']:
            # The value of the dictionary is the "is_obselete" field for that
            # reaction
            potential_rxns[rxn['id']] = rxn['is_obsolete']
    # If there is only one match, use it
    if len(potential_rxns) == 1:
        print("Success: " + bigg_id)
        modelseed_rxn = list(potential_rxns.keys())[0]
    # If there are multiple matches, we may still be able to pick one
    elif len(potential_rxns) > 1:
        # See if any of the matches are not obsolete
        non_obsolete = [x for x in potential_rxns if not potential_rxns[x]]
        # If there is only one non-obsolete match, use it
        if len(non_obsolete) == 1:
            modelseed_rxn = non_obsolete[0]
        # If there are multiple non-obsolete matches, give a warning and skip it
        elif len(non_obsolete) > 1:
            print("WARNING: Multiple non-obselete matches for " + bigg_id)
            continue
    # Change the reaction name (bigg_id) to the ModelSeed ID
    rxn_info['bigg_id'] = modelseed_rxn
    # Change the metabolite IDs in the reaction
    for metabolite in rxn_info['metabolites']:
        # If the metabolite is in the dictionary, change it
        if metabolite['bigg_id'] in bigg_to_modelseed.keys():
            metabolite['bigg_id'] = bigg_to_modelseed[metabolite['bigg_id']]
        # If it isn't, give a warning
        else:
            print("WARNING: " + metabolite['bigg_id'] + " not in dictionary")

# Save the new map
with open('escher/iJO1366.Central metabolism.modelseed.json', 'w') as f:
    json.dump(map, f, indent=2)