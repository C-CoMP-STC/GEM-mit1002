import os

import cobra
import pandas as pd


def search_metabolites_in_model():
    """
    Loads exometabolite data, maps metabolite names to IDs, and searches for them
    in a metabolic model, outputting a sorted table by max concentration.
    """
    # Get the directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Define file paths
    data_file = os.path.join(script_dir, "ProDiel_quant_20260211.csv")
    map_file = os.path.join(script_dir, "metabolite_id_map.csv")
    model_file = os.path.join(
        script_dir, "../../model.xml"
    )  # Go up two directories to find model.xml
    output_file = os.path.join(script_dir, "ex_metabolite_search_results.csv")

    # 1. Load "EX" metabolites
    df_data = pd.read_csv(data_file)
    df_ex = df_data[df_data["INorEX"] == "EX"].copy()

    # 2. Get metabolite IDs
    df_map = pd.read_csv(map_file)

    # Merge dataframes to get IDs
    df_merged = pd.merge(df_ex, df_map, left_on="CleanName", right_on="name")

    # 3. Search in model
    model = cobra.io.read_sbml_model(model_file)
    model_metabolite_ids = {met.id for met in model.metabolites}

    results = []

    # Group by metabolite to process each one
    for name, group in df_merged.groupby("CleanName"):
        # 4. Calculate maximum concentration
        max_conc = group["nM"].max()

        # Assuming one ID per name for simplicity, take the first one
        met_id = group["id"].iloc[0]

        # Check for metabolite in cytoplasm and extracellular compartments
        found_c = f"{met_id}_c0" in model_metabolite_ids
        found_e = f"{met_id}_e0" in model_metabolite_ids

        results.append(
            {
                "MetaboliteName": name,
                "MetaboliteID": met_id,
                "MaxConcentration_nM": max_conc,
                "Found_in_Cytoplasm_c0": found_c,
                "Found_in_Extracellular_e0": found_e,
            }
        )

    # 5. Output a sorted table
    df_results = pd.DataFrame(results)
    df_results_sorted = df_results.sort_values(
        by="MaxConcentration_nM", ascending=False
    )

    # Save the results to a CSV file
    df_results_sorted.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    search_metabolites_in_model()
