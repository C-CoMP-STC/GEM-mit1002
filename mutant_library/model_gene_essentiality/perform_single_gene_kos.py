import os
import warnings

import cobra

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(os.path.dirname(FILE_DIR))
OUT_DIR = os.path.join(FILE_DIR, "results")
os.makedirs(OUT_DIR, exist_ok=True)


def main():
    """Loads model, performs gene deletions, and saves results."""
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(REPO_DIR, "model.xml"))

    # Set the model's media (e.g., minimal media)
    glc_medium = {
        "EX_cpd00027_e0": 10,  # D-Glucose_e0
        "EX_cpd00007_e0": 1000,  # O2_e0
        # Remaining minimal media components
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
        "EX_cpd00075_e0": 1000,  # Nitrite_e0
        "EX_cpd00001_e0": 1000,  # H2O_e0
        "EX_cpd00635_e0": 1000,  # Cbl_e0
        "EX_cpd00034_e0": 1000,  # Zn2+_e0
        "EX_cpd00149_e0": 1000,  # Co2+_e0
    }
    model.medium = clean_media(model, glc_medium)

    # Perform single gene knockouts
    print("Performing single gene knockouts...")
    ko_results = cobra.flux_analysis.single_gene_deletion(model)
    print("Finished knockouts.")

    # Unpack the "ids" results from sets of genes to just the gene IDs
    ko_results["ids"] = ko_results["ids"].apply(
        lambda x: next(iter(x)) if isinstance(x, set) and len(x) == 1 else x
    )

    # Save the results to a CSV file
    output_file = os.path.join(OUT_DIR, "single_gene_ko_results.csv")
    ko_results.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")


# Helper function for setting the media regardless if the exchange reaction is
# present in the model
# TODO: Move this to a helper file
def clean_media(model, media):
    """clean_media
    Removes exchange reactions from the media that are not present in the model

    Parameters
    ----------
    model : cobra.Model
        The model to set the media for.
    media : dict
        A dictionary where the keys are the exchange reactions for the metabolites
        in the media, and the values are the lower bound for the exchange reaction.

    Returns
    -------
    dict
        A dictionary where the keys are the exchange reactions for the metabolites
        in the media, and the values are the lower bound for the exchange reaction
    """
    # Make an empty dictionary for the media
    clean_medium = {}
    # Loop through the media and set the exchange reactions that are present
    for ex_rxn, lb in media.items():
        if ex_rxn in [r.id for r in model.reactions]:
            clean_medium[ex_rxn] = lb
        else:
            warnings.warn(
                "Model does not have the exchange reaction "
                + ex_rxn
                + ", so it was not set in the media."
            )

    # Return the clean medium
    return clean_medium


if __name__ == "__main__":
    main()
