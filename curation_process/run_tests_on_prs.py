import os
import pickle
import subprocess
import warnings

import cobra
import pandas as pd

# Remote name
REMOTE = "origin"

# FILE PATHS
FILE_PATH = os.path.dirname(os.path.abspath(__file__))
REPO_PATH = os.path.dirname(FILE_PATH)
TESTFILE_DIR = os.path.join(REPO_PATH, "test", "test_files")

# Load the media definitions
with open(os.path.join(TESTFILE_DIR, "media", "media_definitions.pkl"), "rb") as f:
    media_definitions = pickle.load(f)
minimal_media = media_definitions["minimal"]


def run_tests_on_prs():
    # Save current branch
    original_branch = subprocess.check_output(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True
    ).strip()

    # List of PRs to evaluate
    # For the initial model construction to growing on everything use PRs 89 - 212
    # To highlight the acetate/leucine/isolecuine fixes use PRs 273-293
    # For initial to after fixing the TCA cycle fluxes use PRs (89-344)
    # but skip merges to the main branch (89, 91, 92, 93, 100, 107, 108, 109, 110, 112, 114, 115, 117, 203, 267, 273)
    pull_requests = list(range(89, 344))
    skip_prs = [89, 91, 92, 93, 100, 107, 108, 109, 110, 112, 114, 115, 117, 203, 267, 273]
    pull_requests = [pr for pr in pull_requests if pr not in skip_prs]

    # Prepare results as a list of dicts
    results_list = []

    for pr in pull_requests:
        branch_name = f"pr-{pr}"
        print(f"\n--- Evaluating PR #{pr} ---")

        try:
            # Delete the branch if it already exists locally
            subprocess.run(
                ["git", "branch", "-D", branch_name],
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

            # Fetch and checkout PR
            subprocess.run(
                ["git", "fetch", REMOTE, f"pull/{pr}/head:{branch_name}"],
                check=True,
            )
            subprocess.run(["git", "checkout", branch_name], check=True)

            # Get commit hash
            commit_hash = subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip()

            # Run test_growth
            results = run_test_growth()

            # Append results to list
            results_list.append(
                {
                    "PR Number": pr,
                    "Commit": commit_hash,
                    "Reactions": results.get("num_reactions"),
                    "Metabolites": results.get("num_metabolites"),
                    "Genes": results.get("num_genes"),
                    "Matches": results.get("matches"),
                    "Total": results.get("total"),
                    "% Match": round(
                        100
                        * results.get("matches", 0)
                        / max(results.get("total", 1), 1),
                        2,
                    ),
                }
            )

        except Exception as e:
            print(f"Error with PR #{pr}: {e}")
            results_list.append(
                {
                    "PR Number": pr,
                    "Commit": "ERROR",
                    "Reactions": "ERROR",
                    "Metabolites": "ERROR",
                    "Genes": "ERROR",
                    "Matches": "ERROR",
                    "Total": "ERROR",
                    "% Match": "ERROR",
                }
            )

    # Return to original branch
    subprocess.run(["git", "checkout", original_branch])

    # Write results to CSV ONCE, after all PRs are processed
    summary_file = os.path.join(FILE_PATH, "growth_match_summary.csv")
    df = pd.DataFrame(results_list)
    df.to_csv(summary_file, index=False)


def run_test_growth():
    model = cobra.io.read_sbml_model(os.path.join(REPO_PATH, "model.xml"))

    # Load the TSV of the growth phenotypes
    growth_phenotypes = pd.read_csv(
        os.path.join(TESTFILE_DIR, "known_growth_phenotypes.tsv"), sep="\t"
    )

    # Count model components
    num_reactions = len(model.reactions)
    num_metabolites = len(model.metabolites)
    num_genes = len(model.genes)

    matches = 0
    total = 0

    for _, row in growth_phenotypes.iterrows():
        # Get the expected growth phenotype
        expected_growth = row["growth"]

        # If the expected growth is not Yes or No (e.g. "Unsure"), skip the row
        if expected_growth not in ["Yes", "No"]:
            warnings.warn(
                f"Expected growth phenotype '{expected_growth}' is not valid. Skipping row."
            )
            continue

        # Set the minimal media and exchange reactions
        minimal_media = media_definitions[row["minimal_media"]].copy()

        # Check if the model has an exchange reaction for the metabolite
        if "EX_" + row["met_id"] + "_e0" in [r.id for r in model.reactions]:
            # If it does, add the exchange reaction to the minimal media used
            minimal_media["EX_" + row["met_id"] + "_e0"] = 1000.0
            # Set the media
            model.medium = clean_media(model, minimal_media)
            # Run the model
            sol = model.optimize()
            # Check if the model grows
            if sol.objective_value > 1e-3:
                # If it does, set to Yes
                pred_growth = "Yes"
            else:
                # If it doesn't, set to No
                pred_growth = "No"
        else:
            # If the exchange reaction is not present, we cannot test growth
            pred_growth = "No"

        if pred_growth == expected_growth:
            matches += 1
        total += 1

    return {
        "matches": matches,
        "total": total,
        "num_reactions": num_reactions,
        "num_metabolites": num_metabolites,
        "num_genes": num_genes,
    }


def clean_media(model: cobra.Model, media: dict) -> dict:
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
    run_tests_on_prs()
    print("Test results saved to growth_match_summary.csv")
    print("You can plot the results using the provided plotting script.")
