import os

import cobra
from gem_utilities.maps import map_ko_ids

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(FILE_DIR, 'kegg_maps')

# If the outpath doesn't exist, create it
os.makedirs(OUT_DIR, exist_ok=True)

# Load the model
model = cobra.io.read_sbml_model(os.path.join(FILE_DIR, 'model.xml'))

# Extract the KO IDs from the model
ko_ids = set()
for reaction in model.reactions:
    for ko_id in reaction.annotation.get('kegg.orthology', []):
        ko_ids.add(ko_id)

# Make and save a single map
# TODO: Add more maps
# TODO: Find all maps to make based on the KO IDs
map_ko_ids("ko00061",
           ko_ids,
           kgml_folder='/Users/helenscott/Documents/PhD/Segre-lab/kegg_data/kgml',
           output_folder=OUT_DIR)
