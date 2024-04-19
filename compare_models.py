import cobra
import pandas as pd

# Load the two models
modelseedpy_model = cobra.io.read_sbml_model("modelseedpy_model.xml")
kbase_model = cobra.io.read_sbml_model("model.xml")

# Make a set of the reactions IDs for any reaction that is in either model
modelseedpy_reactions = set([reaction.id for reaction in modelseedpy_model.reactions])
kbase_reactions = set([reaction.id for reaction in kbase_model.reactions])
union_reactions = modelseedpy_reactions.union(kbase_reactions)

# Make a pandas dataframe of the union reactions with the reacion ID,
# if the reaction is present in each model, the reaction name,
# and the reaction equation
reaction_data = []
for reaction_id in union_reactions:
    row = {"ID": reaction_id}
    # Check if the reaction is in each model
    if reaction_id in kbase_reactions:
        row["KBase Model"] = "Y"
    else:
        row["KBase Model"] = "N"
    if reaction_id in modelseedpy_reactions:
        row["ModelSEEDpy Model"] = "Y"
    else:
        row["ModelSEEDpy Model"] = "N"
    # Get the reaction name and equation
    if reaction_id in modelseedpy_reactions:
        reaction = modelseedpy_model.reactions.get_by_id(reaction_id)
        row["Name"] = reaction.name
        row["Equation"] = reaction.reaction
    else:
        reaction = kbase_model.reactions.get_by_id(reaction_id)
        row["Name"] = reaction.name
        row["Equation"] = reaction.reaction
    reaction_data.append(row)

# Sort the list, so that the reactions only in the KBase model are
# at the top, then the reactions only in the ModelSEEDpy model, then the
# reactions in both models
reaction_df = pd.DataFrame(reaction_data)
reaction_df = reaction_df.sort_values(
    by=["ModelSEEDpy Model", "KBase Model", "ID"], ascending=[True, True, True]
)

# Save the dataframe to a CSV file
reaction_df.to_csv("model_comparison.csv", index=False)
