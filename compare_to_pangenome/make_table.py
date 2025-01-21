import pandas as pd
import cobra


# define a function that returns the MetaNetX ID for a different database's ID
def get_mnx_id(id, source):
    # Combine the source and ID to match the format in the xref file
    combo_id = source + ':' + id
    # Check if there is a match in the xref file, if not, do nothing
    if combo_id not in mnx_xref['source'].values:
        return None
    # If there is more than one match, return all of them
    return mnx_xref[mnx_xref['source'] == combo_id]['ID'].values

def get_seed_id(mnx_id):
    # Check if there is a match in the xref file, if not, do nothing
    # This should always be true, since I just got the ID from this database
    if mnx_id not in mnx_xref['ID'].values:
        return None
    # Get the rows that match the MetaNetX ID
    other_ids = mnx_xref[mnx_xref['ID'] == mnx_id]
    # Get just the IDs that start with "seed.reaction:" for the source
    seed_xrefs = other_ids[other_ids['source'].str.startswith('seed.reaction:')]
    # Remove the "seed.reaction:" from the source and return all of the IDs
    return seed_xrefs['source'].str.replace('seed.reaction:', '').values


# Define a function that searches the values of a dictionary (which are both
# strings and lists) for a string, and return the key(s) that have it as a
# value
def search_dict(dictionary, search):
    keys = []
    for key, value in dictionary.items():
        if search in value:
            keys.append(key)
    return keys


# Read in Michelle's pangenome table
db = pd.read_csv('Pangenome from Michelle/Database_MIT1002GeneCalls.csv')

# Load in my KBase model
model = cobra.io.read_sbml_model("model.xml")

# Read in the MetaNetX reactions xref spreadsheet
mnx_xref = pd.read_csv('reac_xref.tsv', sep='\t', comment='#', header=None)
mnx_xref.columns = ['source', 'ID', 'description']

# Add a column for the MetaNetX reaction ID
db['MetaNetX Reaction ID'] = None

# Get the MetaNetX ID for each reaction in the database
for index, row in db.iterrows():
    if pd.isnull(row['KEGG.Reaction']):
        continue
    mnx_id = get_mnx_id(row['KEGG.Reaction'], 'kegg.reaction')
    # If there is no MetaNetX ID, do nothing
    if not mnx_id:
        continue
    # If there is only one MetaNetX ID, add it to the annotation
    if len(mnx_id) == 1:
        db.at[index, 'MetaNetX Reaction ID'] = mnx_id[0]
    # If there are multiple MetaNetX IDs, add them all as a list
    if len(mnx_id) > 1:
        db.at[index, 'MetaNetX Reaction ID'] = [id for id in mnx_id]

# Get the Seed ID for each reaction in the database
db['ModelSEED ID'] = None
for index, row in db.iterrows():
    if pd.isnull(row['MetaNetX Reaction ID']):
        continue
    seed_id = get_seed_id(row['MetaNetX Reaction ID'])
    # If there is no ModelSEED ID, do nothing
    if len(seed_id) == 0:
        continue
    # Otherwise, add the whole list to the database
    # Have to make a python list, rather than a numpy array, to save it as a list in the csv
    db.at[index, 'ModelSEED ID'] = list(seed_id)

# Save the database with the new columns
db.to_csv('Pangenome from Michelle/database_w_MNX_SEED.csv', index=False)

# # Load in the cleaned up database
# db = pd.read_csv('Pangenome from Michelle/database_w_MNX.csv')

# # Add a column for presence in the KBase model
# db['In KBase Model'] = 0
# db['ModelSEED ID'] = None
# db['Gene Call in KBase Model'] = None

# # Make a dictionary of the KBase Model's reaction where the key is the reaction ID
# # and the value is the MetaNetX ID from the annotation
# model_reactions = {}
# for reaction in model.reactions:
#     if 'metanetx.reaction' in reaction.annotation:
#         model_reactions[reaction.id] = reaction.annotation['metanetx.reaction']

# # Check if each reaction in the database is in the dictionary of the KBase model's reactions
# # and if it is, mark it as present in the KBase model, and save the reaction's
# # ID and gene call
# for index, row in db.iterrows():
#     if pd.isnull(row['MetaNetX Reaction ID']):
#         continue
#     seed_ids = search_dict(model_reactions, row['MetaNetX Reaction ID'])
#     if len(seed_ids) > 1:
#         print('HELEN HELP')
#     if not seed_ids:
#         continue
#     db.at[index, 'In KBase Model'] = 1
#     db.at[index, 'ModelSEED ID'] = seed_ids
#     db.at[index, 'Gene Call in KBase Model'] = [model.reactions.get_by_id(seed_id).gene_reaction_rule for seed_id in seed_ids]

# # Save the database with the new columns
# db.to_csv('Pangenome from Michelle/database_w_KBase.csv', index=False)

# Load in the database with the KBase model information
db = pd.read_csv('Pangenome from Michelle/database_w_KBase.csv')

# Do the same comparison with the ModelSEEDpy model
# Load in the ModelSEEDpy model
modelseedpy_model = cobra.io.read_sbml_model("modelseedpy_model.xml")

# Add a column for presence in the KBase model
db['In ModelSEEDpy Model'] = 0
db['Gene Call in ModelSEEDpy Model'] = None

# Make a dictionary of the ModelSEEDpy Model's reaction where the key is the reaction ID
# and the value is the MetaNetX ID from the annotation
msp_model_reactions = {}
for reaction in modelseedpy_model.reactions:
    if 'metanetx.reaction' in reaction.annotation:
        msp_model_reactions[reaction.id] = reaction.annotation['metanetx.reaction']

# Check if each reaction in the database is in the dictionary of the KBase model's reactions
# and if it is, mark it as present in the KBase model, and save the reaction's
# ID and gene call
for index, row in db.iterrows():
    if pd.isnull(row['MetaNetX Reaction ID']):
        continue
    seed_ids = search_dict(msp_model_reactions, row['MetaNetX Reaction ID'])
    if len(seed_ids) > 1:
        print('HELEN HELP')
    if not seed_ids:
        continue
    db.at[index, 'In ModelSEEDpy Model'] = 1
    db.at[index, 'Gene Call in ModelSEEDpy Model'] = [modelseedpy_model.reactions.get_by_id(seed_id).gene_reaction_rule for seed_id in seed_ids]

# Save the database with the new columns
db.to_csv('Pangenome from Michelle/database_w_MSP.csv', index=False)