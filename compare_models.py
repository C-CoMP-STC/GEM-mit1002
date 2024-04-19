import cobra

# Load the two models
modelseedpy_model = cobra.io.read_sbml_model("modelseedpy_model.xml")
kbase_model = cobra.io.read_sbml_model("model.xml")

# Make a list of the reactions that the two models do not share
modelseedpy_unique_reactions = []
kbase_unique_reactions = []
for reaction in modelseedpy_model.reactions:
    if reaction.id not in kbase_model.reactions:
        modelseedpy_unique_reactions.append(reaction)
for reaction in kbase_model.reactions:
    if reaction.id not in modelseedpy_model.reactions:
        kbase_unique_reactions.append(reaction)

# Save the unique reactions to a file
with open("modelseedpy_unique_reactions.txt", "w") as f:
    for reaction in modelseedpy_unique_reactions:
        f.write(reaction.id + "\n")
with open("kbase_unique_reactions.txt", "w") as f:
    for reaction in kbase_unique_reactions:
        f.write(reaction.id + "\n")
