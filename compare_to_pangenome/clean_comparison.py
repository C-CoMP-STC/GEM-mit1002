import pandas as pd

# Load the database with the comparison to the ModelSEEDpy model
db = pd.read_csv("Pangenome from Michelle/database_w_MSP.csv")

# Filter the database to only have rows that are present in the pangenome
# (1 in the "Inferred Presenece" column) OR are present in the ModelSEEDpy model
# (1 in the "In ModelSEEDpy Model" column)
db_present_reactions = db[
    (db["Inferred Presence"] == 1) | (db["In ModelSEEDpy Model"] == 1)
]

# Make a new empty dataframe
clean_db = pd.DataFrame(
    columns=[
        "MetaNetX Reaction ID",
        "KEGG ID",
        "SEED ID",
        "Inferred in Pangenome",
        "Pangenome Gene Call",
        "Present in ModelSEEDpy Model",
        "ModelSEEDpy Gene Call",
    ]
)

# Get a list of the unique MetaNetX Reaction IDs
unique_metanetx_ids = db_present_reactions["MetaNetX Reaction ID"].unique()

# Loop through the unique MetaNetX Reaction IDs
for metanetx_id in unique_metanetx_ids:
    # Make a dictionary for the new row's information
    new_row = {"MetaNetX Reaction ID": metanetx_id}

    # Get the rows that have the current MetaNetX Reaction ID
    rows = db_present_reactions[
        db_present_reactions["MetaNetX Reaction ID"] == metanetx_id
    ]

    # Skip if there are no rows (happens with metanetx id is nan)
    if len(rows) == 0:
        continue

    # Get the KEGG ID, SEED ID, and presence from the rows (check that they are all the same)
    kegg_ids = rows["KEGG.Reaction"].unique()
    if len(kegg_ids) != 1:
        new_row["KEGG ID"] = 'Inconsistent'
    elif len(kegg_ids) == 1:
        new_row["KEGG ID"] = kegg_ids[0]
    seed_ids = rows["ModelSEED ID"].unique()
    if len(seed_ids) != 1:
        new_row["SEED ID"] = 'Inconsistent'
    elif len(seed_ids) == 1:
        new_row["SEED ID"] = seed_ids[0]
    inferred_presence = rows["Inferred Presence"].unique()
    if len(inferred_presence) != 1:
        new_row["Inferred in Pangenome"] = 'Inconsistent'
    elif len(inferred_presence) == 1:
        new_row["Inferred in Pangenome"] = inferred_presence[0]
    modelseedpy_presence = rows["In ModelSEEDpy Model"].unique()
    if len(modelseedpy_presence) != 1:
        new_row["Present in ModelSEEDpy Model"] = 'Inconsistent'
    elif len(modelseedpy_presence) == 1:
        new_row["Present in ModelSEEDpy Model"] = modelseedpy_presence[0]
    modelseed_gene_calls = rows["Gene Call in ModelSEEDpy Model"].unique()
    if len(modelseed_gene_calls) != 1:
        new_row["ModelSEEDpy Gene Call"] = 'Inconsistent'
    elif len(modelseed_gene_calls) == 1:
        new_row["ModelSEEDpy Gene Call"] = modelseed_gene_calls[0]

    # Start a list for the gene-reaction rule
    gene_reaction_rule = []
    # Get a list of unique "Enzyme" values
    unique_enzymes = rows["Enzyme"].unique()
    # Loop thorugh the unique "Enzyme" values
    for enzyme in unique_enzymes:
        # Get the rows that have the current "Enzyme" value
        enzyme_rows = rows[rows["Enzyme"] == enzyme]
        # If there is only one row, get the gene calls
        if len(enzyme_rows) == 1:
            gene_calls = enzyme_rows["Gene.Call.ID"].to_string(index=False)
        # If there are multiple rows, concatenate the gene calls together with 'and' and put it all in parentheses
        elif len(enzyme_rows) > 1:
            gene_calls = "("
            for index, row in enzyme_rows.iterrows():
                if index == 0:
                    gene_calls += str(row["Gene.Call.ID"])
                else:
                    gene_calls += " and " + str(row["Gene.Call.ID"])
            gene_calls += ")"
        # Add that enzyme's gene calls to the gene-reaction rule string
        gene_reaction_rule.append(gene_calls)
    # Join the gene calls together with ' or ' and put it all in parentheses
    gene_reaction_rule_string = "['" + " or ".join(gene_reaction_rule) + "']"
    # Add the gene-reaction rule string to the new row
    new_row["Pangenome Gene Call"] = gene_reaction_rule_string
    clean_db = clean_db.append(new_row, ignore_index=True)

# Save the cleaned database
clean_db.to_csv("Pangenome from Michelle/clean_compariosn.csv", index=False)
