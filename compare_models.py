import cobra

# Load the two models
modelseedpy_model = cobra.io.read_sbml_model("modelseedpy_model.xml")
kbase_model = cobra.io.read_sbml_model("model.xml")

# Make a set of the reactions IDs for any reaction that is in either model
modelseedpy_reactions = set([reaction.id for reaction in modelseedpy_model.reactions])
kbase_reactions = set([reaction.id for reaction in kbase_model.reactions])
union_reactions = modelseedpy_reactions.union(kbase_reactions)

# Make a table of the union reactions with the reacion ID,
# if the reaction is present in each model, the reaction name,
# and the reaction equation
with open("union_reactions.tsv", "w") as f:
    f.write("ID\tKBase Model\tModelSEEDpy Model\tName\tEquation\n")
    for reaction_id in union_reactions:
        f.write(f"{reaction_id}\t")
        # Check if the reaction is in each model
        if reaction_id in kbase_reactions:
            f.write("Y\t")
        else:
            f.write("\t")
        if reaction_id in modelseedpy_reactions:
            f.write("Y\t")
        else:
            f.write("\t")
        # Get the reaction name and equation
        if reaction_id in modelseedpy_reactions:
            reaction = modelseedpy_model.reactions.get_by_id(reaction_id)
            f.write(f"{reaction.id}\t{reaction.name}\t{reaction.reaction}\n")
        else:
            reaction = kbase_model.reactions.get_by_id(reaction_id)
            f.write(f"{reaction.id}\t{reaction.name}\t{reaction.reaction}\n")
