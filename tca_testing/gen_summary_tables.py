#!/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/GEM-mit1002/.venv/bin/python
import ast

import pandas as pd

# Load results
orig_df = pd.read_csv("results/amac_Original_growth_fluxes.csv")
strict_df = pd.read_csv("results/amac_Strict_ATP_Production_growth_fluxes.csv")
rescue_df = pd.read_csv("results/amac_Strict_ATP_Production_With_dUTP_Nucleotidohydrolase_growth_fluxes.csv")

# Table of growth rate and ATP-consuming reactions with negative flux
# Convert the 'fluxes' column from string representation of dict to actual dict
orig_df["fluxes"] = orig_df["fluxes"].apply(ast.literal_eval)

atp_consuming_reactions = [
     "rxn00077_c0",
     "rxn00104_c0",
     "rxn00239_c0",
     "rxn00364_c0",
     "rxn00379_c0",
     "rxn01219_c0",
     "rxn01509_c0",
     "rxn01517_c0",
     "rxn02314_c0",
     "rxn08762_c0",
     "rxn15121_c0",
]

# Create a new column that lists the specific reactions with negative flux
orig_df["negative_atp_flux_rxns"] = orig_df[
    "fluxes"
].apply(
    lambda flux_dict: [
        rxn
        for rxn in atp_consuming_reactions
        if flux_dict.get(rxn, 0) < -1e-9  # Use a small tolerance
    ]
)

# Filter for rows where the list of negative flux reactions is not empty
rows_with_negative_flux = orig_df[
    orig_df["negative_atp_flux_rxns"].str.len() > 0
]

# Select the most relevant columns for the output
display_cols = [
    "Model",
    "C_source",
    "N_source",
    "growth_rate",
    "negative_atp_flux_rxns",
]

# Print the resulting dataframe in markdown format
print(rows_with_negative_flux[display_cols].to_markdown(index=False))

# Table comparing growth rates across conditions
comparison_df = pd.merge(
    orig_df[["Model", "C_source", "N_source", "growth_rate"]],
    strict_df[["C_source", "N_source", "growth_rate"]],
    on=["C_source", "N_source"],
    suffixes=("_Original", "_Strict_ATP_Production"),
)
comparison_df = pd.merge(
    comparison_df,
    rescue_df[["C_source", "N_source", "growth_rate"]],
    on=["C_source", "N_source"],
)
comparison_df = comparison_df.rename(
    columns={"growth_rate": "growth_rate_With_dUTP_Nucleotidohydrolase"}
)
# Print the comparison dataframe in markdown format
print(comparison_df.to_markdown(index=False))