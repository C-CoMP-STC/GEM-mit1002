import cobra
import json
import pandas as pd

# Load the model from the SBML file
model = cobra.io.read_sbml_model('model.xml')

# Export the model to JSON
cobra.io.save_json_model(model, 'model.json')

# Convert to excel file
# Load the json model 
with open('model.json') as f:
    model_json = json.load(f)

met_df = pd.DataFrame(model_json['metabolites'])
rxn_df = pd.DataFrame(model_json['reactions'])
gene_df = pd.DataFrame(model_json['genes'])

# Save to excel
with pd.ExcelWriter('model.xlsx') as writer:
    met_df.to_excel(writer, sheet_name='metabolites', index=False)
    rxn_df.to_excel(writer, sheet_name='reactions', index=False)
    gene_df.to_excel(writer, sheet_name='genes', index=False)
