import cobra
import pickle as pkl

# Load the model
model = cobra.io.read_sbml_model("model.xml")

# Load the media
with open("test/test_files/media/media_definitions.pkl", "rb") as f:
    media = pkl.load(f)

# Add glucose and remove cobalamin (vitamin with C) from the media
default_medium = media['minimal'].copy()
default_medium['EX_cpd00027_e0'] = 10
default_medium['EX_cpd00635_e0'] = 0

# Set the medium
model.medium = default_medium

# Save the model
cobra.io.write_sbml_model(model, "model.xml")