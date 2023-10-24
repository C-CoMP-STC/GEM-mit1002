import cobra

# Load the model from the SBML file
model = cobra.io.read_sbml_model('model.xml')

# Export the model to JSON
cobra.io.save_json_model(model, 'model.json')
