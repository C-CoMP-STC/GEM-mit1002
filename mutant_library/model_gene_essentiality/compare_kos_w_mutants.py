import os

import pandas as pd

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(FILE_DIR, "results")
MUTANT_LIB_DIR = os.path.dirname(FILE_DIR)

# Define the file paths for your input files
lookup_file = os.path.join(MUTANT_LIB_DIR, "MIT1002_mutant-library-gene-lookup.csv")
mutant_lib_file = os.path.join(MUTANT_LIB_DIR, "MIT2001-735-unique-colonies.csv")
results_file = os.path.join(OUT_DIR, "single_gene_ko_results.csv")
output_file = os.path.join(OUT_DIR, "merged_ko_results_with_mut_lib_id.csv")

# Define the column names for the merge keys for merging the gene IDs
# This is the column from the lookup table
lookup_key_column = "locus_tag_from_model"
# This is the column from your results file that matches the lookup key
results_key_column = "ids"

# Define the column names for merging the gene IDs with the mutant library data
mutant_lib_key_column = "locus_tag"

# --- Main Script ---
try:
    # Load the CSV files into pandas DataFrames
    print(f"Reading lookup data from '{lookup_file}'...")
    lookup_df = pd.read_csv(lookup_file)

    print(f"Reading single gene KO results from '{results_file}'...")
    results_df = pd.read_csv(results_file)

    print(f"Reading mutant library data from '{mutant_lib_file}'...")
    mutant_lib_df = pd.read_csv(mutant_lib_file)

    # Prepare the lookup table for merging
    # We only need the two columns for the mapping.
    # This also handles cases where one 'locus_tag_from_model' might have multiple 'locus_tag_from_mut_lib'
    # by keeping the first one found.
    mapping_df = lookup_df[
        [lookup_key_column, "locus_tag_from_mut_lib"]
    ].drop_duplicates(subset=[lookup_key_column])

    # Add a new column to the results dataframe with boolean values indicating if the gene is essential in the model
    # A gene is considered essential if the growth rate is less that 1E-3
    results_df["essential_in_model"] = results_df["growth"] < 1e-3

    # Merge the two DataFrames with gene IDs
    # We use a 'left' merge to ensure all rows from the results_df are kept.
    print(
        f"Merging tables on '{results_key_column}' (from results) and '{lookup_key_column}' (from lookup)..."
    )
    merged_df = pd.merge(
        results_df,
        mapping_df,
        left_on=results_key_column,
        right_on=lookup_key_column,
        how="outer",
    )

    # Drop the redundant key column from the lookup table if you don't need it
    if (
        lookup_key_column in merged_df.columns
        and lookup_key_column != results_key_column
    ):
        merged_df = merged_df.drop(columns=[lookup_key_column])

    # Merge with the mutant library data
    print(
        f"Merging with mutant library data on '{mutant_lib_key_column}' (from mutant lib) and 'locus_tag_from_mut_lib' (from lookup)..."
    )
    merged_df = pd.merge(
        merged_df,
        mutant_lib_df,
        left_on="locus_tag_from_mut_lib",
        right_on=mutant_lib_key_column,
        how="left",
    )

    # Save the new DataFrame to a CSV file
    merged_df.to_csv(output_file, index=False)

    print("-" * 20)
    print("Merge complete!")
    print(f"Output saved to '{output_file}'")

except FileNotFoundError as e:
    print(f"Error: {e}. Please make sure both CSV files are in the correct location.")
except KeyError as e:
    print(f"Error: A column name was not found: {e}.")
    print(
        f"Please check that '{results_key_column}' is the correct column name in '{results_file}'."
    )
