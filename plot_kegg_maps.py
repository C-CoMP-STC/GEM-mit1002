from gem_utilities.maps import map_ko_ids

# TODO: Replace by extracting all of the KOs from the model
ko_ids = ["K00665", "K00668"]

# Make and save a single map
# TODO: Change the KGML folder
# TODO: Change the output folder
# TODO: Add more maps
# TODO: Find all maps to make based on the KO IDs
map_ko_ids("ko00061", ko_ids, kgml_folder='.')