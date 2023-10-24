import cobra
import pandas as pd

# Load the TSV of the growth phenotypes
growth_phenotypes = pd.read_csv("known_growth_phenotypes.tsv", sep="\t")

# Rename the experimental column
growth_phenotypes = growth_phenotypes.rename(columns={"growth": "exp_growth"})

# Load the model
model = cobra.io.read_sbml_model("model.xml")

# Make a minimal media (not including any carbon sources)
no_c_media = {
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00007_e0": 20,  # O2_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00075_e0": 1000,  # Nitrite_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# Loop through the growth phenotypes and add the predicted phenotype to a list
pred_growth = []
for index, row in growth_phenotypes.iterrows():
    # Check if the model has an exchange reaction for the metabolite
    if "EX_" + row["met_id"] + "_e0" in [r.id for r in model.reactions]:
        # If it does, add the exchange reaction to the minimal media
        minimal_media = no_c_media.copy()
        minimal_media["EX_" + row["met_id"] + "_e0"] = 1000.0
        # Set the media
        model.medium = minimal_media
        # Run the model
        sol = model.optimize()
        # Check if the model grows
        if sol.objective_value > 0:
            # If it does, add 'Y' to the list
            pred_growth.append("Yes")
        else:
            # If it doesn't, add 'N' to the list
            pred_growth.append("No")
    else:
        # If it doesn't have the exchange reaction, add None to the list
        pred_growth.append("No Exchange")

# Add the list as a new column in the dataframe
growth_phenotypes["pred_growth"] = pred_growth

# Save the dataframe as a TSV
growth_phenotypes.to_csv("known_growth_phenotypes_w_pred.tsv", sep="\t", index=False)
