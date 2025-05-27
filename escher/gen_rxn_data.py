import os

import cobra
import pickle

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(FILE_DIR)

# Load the media defintions from the pickle file
with open(os.path.join(REPO_DIR, "test", "test_files", "media", "media_definitions.pkl"), "rb") as f:
    media_definitions = pickle.load(f)

# Load the model
model = cobra.io.read_sbml_model(os.path.join(REPO_DIR, "model.xml"))

####################################
# Glucose as sole carbon source
####################################
# Set the medium to minimal medium
model.medium = media_definitions["minimal_glucose"]
# Run the simulation
sol = model.optimize()
# Save the results as a dictionary with the reaction ID as the key and the
# flux value as the value
glc_rxn_data = {}
