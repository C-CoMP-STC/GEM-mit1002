import json
import os
import pickle
import subprocess
import warnings

import cobra
import pandas as pd
from gem_utilities import media

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


def run_tests_on_prs(pr_start=89, pr_end=None):
    # Save current branch
    original_branch = subprocess.check_output(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True
    ).strip()

    # Get a list of PRs merged into the dev branch (could change to main, by
    # providing target_branch="main")
    all_prs = get_prs_by_target()
    # Filter PRs to those that changed the model.xml file
    model_prs = [pr for pr in all_prs if is_model_changed_in_pr(pr)]
    # Filter PRs for a specific range of interest
    # For the initial model construction to growing on everything use PRs 89 - 212
    # To highlight the acetate/leucine/isolecuine fixes use PRs 273-293
    # For initial to after fixing the TCA cycle fluxes use PRs 89-344
    if pr_end is None:
        pr_end = max(model_prs)
    pull_requests = [pr for pr in model_prs if pr_start <= pr <= pr_end]

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


def run_test_growth(unbounded_flux_limit: int = 1000, biomass_rxn_id="bio1_biomass"):
    """
    On the current branch, this function re-runs the growth with pFBA, while
    holding the amount carbon constant across sources, and saves how many
    media conditions support growth, and how many reactions across all
    simulations have a flux above an arbitrary threshold.

    Parameters
    ----------
    unbounded_flux_limit : int, optional
        The value at which a flux is considered problematically large, by
        default 1000
    biomass_rxn_id : str, optional
        Reaction ID for the biomass reaction, by default "bio1_biomass"

    Returns
    -------
    dict
       Dictionary of the results for the model, including the number of
       reactions, metabolites, and genes, the number of matches between
       predicted and expected growth phenotypes, the total number of phenotypes
       tested, and the number of unique reactions with a flux above the
       unbounded_flux_limit across all simulations.
    """
    model = cobra.io.read_sbml_model(os.path.join(REPO_PATH, "model.xml"))

    # Load the TSV of the growth phenotypes
    growth_phenotypes = pd.read_csv(
        os.path.join(TESTFILE_DIR, "known_growth_phenotypes.tsv"), sep="\t"
    )

    # Count model components
    num_reactions = len(model.reactions)
    num_metabolites = len(model.metabolites)
    num_genes = len(model.genes)

    # Start counters
    matches = 0
    total = 0
    unique_rxns_with_unbounded_flux = set()

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

        # Does this handle multiple meatbolites correctly?
        # Check if the model has an exchange reaction for the metabolite
        if "EX_" + row["met_id"] + "_e0" in [r.id for r in model.reactions]:
            # Get the metabolite object
            met = model.metabolites.get_by_id(row["met_id"] + "_e0")
            # Get the number of carbon atoms in the metabolite
            n_carbons = met.elements.get("C", 0)
            # If the metabolite has carbons, set the lower bound to be 60/n_carbons (equivalent to 10 for glucose)
            # If the metabolite does not have carbon, set an unlimited amount
            if n_carbons == 0:
                # If it does, add the exchange reaction to the minimal media used
                minimal_media["EX_" + row["met_id"] + "_e0"] = 1000.0
            else:
                minimal_media["EX_" + row["met_id"] + "_e0"] = 60 / n_carbons
            # Set the media
            model.medium = media.clean_media(model, minimal_media)
            # Run pFBA on the model
            sol = cobra.flux_analysis.pfba(model)
            # Check if the model grows
            if sol.fluxes[biomass_rxn_id] > 1e-3:
                # If it does, set to Yes
                pred_growth = "Yes"
                # Get the number of reactions in the solution with a flux above the unbounded_flux_limit
                rxns_with_unbounded_flux = [
                    r
                    for r, flux in sol.fluxes.items()
                    if abs(flux) > unbounded_flux_limit
                ]
                # Add the unique reactions with unbounded fluxes to the set
                unique_rxns_with_unbounded_flux.update(rxns_with_unbounded_flux)
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
        "num_reactions": num_reactions,
        "num_metabolites": num_metabolites,
        "num_genes": num_genes,
        "matches": matches,
        "total": total,
        "num_unbounded_rxns": len(unique_rxns_with_unbounded_flux),
    }


def get_prs_by_target(target_branch="dev"):
    # This requires the GitHub CLI 'gh' to be installed and authenticated
    cmd = [
        "gh", "pr", "list",
        "--base", target_branch,
        "--state", "merged",
        "--limit", "1000",
        "--json", "number"
    ]
    output = subprocess.check_output(cmd, text=True)
    return [pr["number"] for pr in json.loads(output)]


def is_model_changed_in_pr(pr_number):
    # Get the list of files changed in the PR
    cmd = ["gh", "pr", "diff", str(pr_number), "--name-only"]
    output = subprocess.check_output(cmd, text=True)
    changed_files = output.splitlines()
    # Check if 'model.xml' is in the list of changed files
    return "model.xml" in changed_files


if __name__ == "__main__":
    # run_tests_on_prs()
    # To debug: Run just the test function on the current branch
    # 10 is way too low, but that way I can check that the counting of unbounded fluxes is working correctly
    results = run_test_growth(unbounded_flux_limit=10)
    print("Test results saved to growth_match_summary.csv")
    print("You can plot the results using the provided plotting script.")
