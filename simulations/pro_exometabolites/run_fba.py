import os

import cobra
import pandas as pd


def run_fba_for_each_timepoint():
    """
    Loads exometabolite data, generates media for each timepoint using Michaelis-Menten kinetics,
    runs pFBA, and saves the results.
    """
    # Get the directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Define file paths
    # TODO: Make this more robust by using relative paths or config files
    data_file = os.path.join(script_dir, "data", "ProDiel_quant_20260211.csv")
    map_file = os.path.join(script_dir, "metabolite_id_map.csv")
    model_file = os.path.join(script_dir, "../../model.xml")
    output_file = os.path.join(script_dir, "results", "fba_results.csv")

    # If the output directory doesn't exist, create it
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Load the model
    model = cobra.io.read_sbml_model(model_file)

    # Load data
    df_data = pd.read_csv(data_file)

    # Filter for "EX" metabolites
    df_ex = df_data[df_data["INorEX"] == "EX"].copy()

    # Load metabolite ID mapping
    df_map = pd.read_csv(map_file)

    # Merge to get metabolite IDs
    df_merged = pd.merge(df_ex, df_map, left_on="CleanName", right_on="name")

    # Michaelis-Menten parameters (default values)
    Vmax = 10  # mmol/gDW/hr
    Km = 1  # mM (assuming nM in data needs conversion)

    results = []

    # Define a minimal medium
    minimal_media = {
        "EX_cpd00007_e0": 20,  # O2_e0
        "EX_cpd00013_e0": 1000,  # Ammonia
        "EX_cpd00058_e0": 1000,  # Cu2+_e0
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
        "EX_cpd00001_e0": 1000,  # H2O_e0
        "EX_cpd00034_e0": 1000,  # Zn2+_e0
        "EX_cpd00149_e0": 1000,  # Co2+_e0
    }

    # Iterate over each timepoint
    for timepoint, group in df_merged.groupby("timepoint"):
        print(f"Processing timepoint: {timepoint}")

        media = minimal_media.copy()

        # Calculate uptake for each metabolite at this timepoint
        for _, row in group.iterrows():
            met_id = row["id"]
            concentration_nM = row["nM"]
            concentration_mM = concentration_nM / 1e6  # Convert nM to mM for Km

            # Michaelis-Menten equation
            uptake_rate = Vmax * concentration_mM / (Km + concentration_mM)

            # Account for "hot spots" within the bulk data by scaling the
            # uptake rate by an arbitrary factor
            hot_spot_factor = 1000
            uptake_rate *= hot_spot_factor

            # Only set the bound if the uptake rate is above a small threshold to avoid numerical issues
            # TODO: Should this be even lower? Maybe 0.001?
            # I did this cause I was getting some very small uptake rates that
            # were causing issues with the solver. But maybe it's better to
            # just set a lower threshold or use a different approach to
            # handle small rates.
            if uptake_rate > 0.01:
                # Find the corresponding exchange reaction
                exchange_reaction_id = f"EX_{met_id}_e0"
                if exchange_reaction_id in model.reactions:
                    media[exchange_reaction_id] = uptake_rate

        # Set the medium for the model
        with model:
            model.medium = media

            # Run pFBA
            try:
                solution = cobra.flux_analysis.pfba(model)
                fluxes = solution.fluxes.to_dict()
            except Exception as e:
                print(f"pFBA failed for timepoint {timepoint}: {e}")
                fluxes = None

        # Store results
        results.append({"time_point": timepoint, "media": media, "fluxes": fluxes})

    # Create and save the results dataframe
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_file, index=False)
    print(f"FBA results saved to {output_file}")


if __name__ == "__main__":
    run_fba_for_each_timepoint()
