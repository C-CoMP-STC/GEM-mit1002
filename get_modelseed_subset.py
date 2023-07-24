import json
import pandas as pd

# Open the model json file
with open('model.json', 'r') as f:
    model = json.load(f)

# Load the ModelSEED reactions tsv as a pandas dataframe
df = pd.read_csv('../../ModelSEEDDatabase/Biochemistry/reactions.tsv',
                 sep='\t')

# Filter the dataframe to only include the reactions in the model
df = df[df['id'].isin([reaction['id'][:-3] for reaction in model['reactions']])]

# Save the filtered dataframe as a tsv
df.to_csv('ALT_modelseed_reactions.tsv', sep='\t', index=False)
