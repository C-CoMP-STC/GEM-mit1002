import pickle as pkl

import cobra

# Load the model
model = cobra.io.read_sbml_model("model.xml")

# Load the media
with open("test/test_files/media/media_definitions.pkl", "rb") as f:
    media = pkl.load(f)

# Add glucose and remove cobalamin (vitamin with C) from the media
default_medium = media["minimal"].copy()
default_medium["EX_cpd00027_e0"] = 10
default_medium["EX_cpd00635_e0"] = 0

# Check if the reaction is there, and if not add it
for r in default_medium.keys():
    if r not in model.reactions:
        # Check if there is an external metabolite with the same ID
        met_id = r.replace(
            "EX_", ""
        )  # Assuming the metabolite ID is the same as the reaction ID without 'EX_'
        try:
            met_e = model.metabolites.get_by_id(met_id)
        except KeyError:
            # If not, add it
            # Get the cytoplasmic metabolite to copy information from
            # This assumes that there is a cytoplasmic version of the metabolite
            # If there weren't, we would have to get the metabolite information from the modelSEED database
            met_c = model.metabolites.get_by_id(met_id.replace("e0", "c0"))
            # Create the metabolite
            met_e = cobra.Metabolite(
                met_id, formula=met_c.formula, charge=met_c.charge, compartment="e0"
            )
            # Add all of the annotations from the cytoplasmic metabolite
            # This assumes that none of the annotations are compartment-specific
            for k, v in met_c.annotation.items():
                met_e.annotation[k] = v
            # Add the metabolite to the model
            model.add_metabolites([met_e])
        # Add the exchange reaction for that metabolite
        model.add_boundary(
            met_e,
            type="exchange",
            reaction_id=r,
            lb=-1 * default_medium[r],
            ub=1000,
        )

# Set the medium
model.medium = default_medium

# Save the model
cobra.io.write_sbml_model(model, "model.xml")
