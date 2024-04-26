import pandas as pd
import cobra

# Load the MetaNetX reactions xref spreadsheet
file_path = 'reac_xref.tsv'
mnx_xref = pd.read_csv(file_path, sep='\t', comment='#', header=None)
mnx_xref.columns = ['source', 'ID', 'description']


# define a function that returns the MetaNetX ID for a different database's ID
def get_mnx_id(id, source):
    # Combine the source and ID to match the format in the xref file
    combo_id = source + ':' + id
    # Check if there is a match in the xref file, if not, do nothing
    if combo_id not in mnx_xref['source'].values:
        return None
    # If there is more than one match, return all of them
    return mnx_xref[mnx_xref['source'] == combo_id]['ID'].values


def add_all_mnx_ids(model, source):
    for reaction in model.reactions:
        if 'metanetx.reaction' not in reaction.annotation:
            # Remove the compartment from the ID
            reaction_id = "_".join(reaction.id.split('_')[0:-1])
            mnx_id = get_mnx_id(reaction_id, source)
            # If there is no MetaNetX ID, do nothing
            if not mnx_id:
                continue
            # If there is only one MetaNetX ID, add it to the annotation
            if len(mnx_id) == 1:
                reaction.annotation['metanetx.reaction'] = mnx_id[0]
            # If there are multiple MetaNetX IDs, add them all as a list
            if len(mnx_id) > 1:
                reaction.annotation['metanetx.reaction'] = [id for id in mnx_id]
    return model


# Load my model files and add the annotation if the MetaNetX ID is currently
# missing
model = cobra.io.read_sbml_model("model.xml")
model = add_all_mnx_ids(model, 'seed.reaction')
cobra.io.write_sbml_model(model, "model.xml")

# Do the same for the ModelSEEDpy model
modelseedpy_model = cobra.io.read_sbml_model("modelseedpy_model.xml")
modelseedpy_model = add_all_mnx_ids(modelseedpy_model, 'seed.reaction')
cobra.io.write_sbml_model(modelseedpy_model, "modelseedpy_model.xml")
