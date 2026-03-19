import json
import os
import pickle
import subprocess
import warnings

import cobra
import pandas as pd
from gem_utilities import media

# FILE PATHS
FILE_PATH = os.path.dirname(os.path.abspath(__file__))
REPO_PATH = os.path.dirname(FILE_PATH)
TESTFILE_DIR = os.path.join(REPO_PATH, "test", "test_files")

# Load the media definitions
with open(os.path.join(TESTFILE_DIR, "media", "media_definitions.pkl"), "rb") as f:
    media_definitions = pickle.load(f)

# Load the TSV of the growth phenotypes
growth_phenotypes = pd.read_csv(
    os.path.join(TESTFILE_DIR, "known_growth_phenotypes.tsv"),
    sep="\t",
    converters={"met_id": lambda x: x.split(",")},
)

def run_tests_on_prs(pr_start=89, pr_end=None):
    # Get a list of PRs merged into the dev branch (could change to main, by
    # providing target_branch="main")
    # By only looking at dev we miss a few PRs that were merged into main,
    # (89-117) but this is easier to automate and still captures the majority
    # of changes to the model
    all_pr_entries = get_prs_by_target()
    # Hardcode PRs to skip because the diff is too large
    # TODO: Instead of hardcoding, change is_model_changed_in_pr to use API
    prs_to_skip = [327]  # Checked manually that the model file was not changed
    # Filter PRs for a specific range of interest
    # For the initial model construction to growing on everything use PRs 89 - 212
    # To highlight the acetate/leucine/isolecuine fixes use PRs 273-293
    # For initial to after fixing the TCA cycle fluxes use PRs 89-344
    # If no upper limit to the range is given, go to the latest version
    if pr_end is None:
        pr_end = max(pr_entry["number"] for pr_entry in all_pr_entries)
    pr_entries_to_check = [
        pr_entry
        for pr_entry in all_pr_entries
        if pr_start <= pr_entry["number"] <= pr_end
        and pr_entry["number"] not in prs_to_skip
    ]
    # Filter PRs to those that changed the model.xml file
    pull_requests = [
        pr_entry
        for pr_entry in pr_entries_to_check
        if is_model_changed_in_pr(pr_entry["number"])
    ]

    # Prepare results as a list of dicts
    results_list = []

    for pr in pull_requests:
        print(f"\n--- Evaluating PR #{pr['number']} ---")

        try:
            # Use gh to get the content of model.xml for a specific PR
            subprocess.run(
                [
                    "gh",
                    "api",
                    f"repos/:owner/:repo/contents/model.xml?ref=pull/{pr['number']}/head",
                    "-H",
                    "Accept: application/vnd.github.v3.raw",
                ],
                stdout=open("temp_model.xml", "w"),
            )

            # Run test_growth
            results = run_test_growth()

            # Append results to list
            results_list.append(
                {
                    "PR Number": pr["number"],
                    "Date Opened": pr["createdAt"],
                    "Date Merged": pr["mergedAt"],
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
                    "Unbounded Flux Reactions": results.get("num_unbounded_rxns"),
                }
            )

        except Exception as e:
            print(f"Error with PR #{pr['number']}: {e}")
            results_list.append(
                {
                    "PR Number": pr["number"],
                    "Date Opened": pr["createdAt"],
                    "Date Merged": pr["mergedAt"],
                    "Reactions": "ERROR",
                    "Metabolites": "ERROR",
                    "Genes": "ERROR",
                    "Matches": "ERROR",
                    "Total": "ERROR",
                    "% Match": "ERROR",
                    "Unbounded Flux Reactions": "ERROR",
                }
            )

    # Delete the temporary model file
    if os.path.exists("temp_model.xml"):
        os.remove("temp_model.xml")

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
    model = cobra.io.read_sbml_model(os.path.join(REPO_PATH, "temp_model.xml"))

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

        # Add the metabolite(s) specified in the row to the media
        for met_id in row["met_id"]:
            # Check that there is an exchange reaction for the metabolite in the model
            if "EX_" + met_id + "_e0" not in [r.id for r in model.reactions]:
                warnings.warn(
                    f"Model does not have an exchange reaction for {met_id}."
                )
                continue
            # Get the metabolite object
            met = model.metabolites.get_by_id(met_id + "_e0")
            # Get the number of carbon atoms in the metabolite
            n_carbons = met.elements.get("C", 0)
            # If the metabolite has carbons, set the lower bound to be 60/n_carbons (equivalent to 10 for glucose)
            # If the metabolite does not have carbon, set an unlimited amount
            if n_carbons == 0:
                minimal_media["EX_" + met_id + "_e0"] = 1000.0
            else:
                minimal_media["EX_" + met_id + "_e0"] = 60 / n_carbons
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
        "gh",
        "pr",
        "list",
        "--base",
        target_branch,
        "--state",
        "merged",
        "--limit",
        "1000",
        "--json",
        "number,createdAt,mergedAt",
    ]
    output = subprocess.check_output(cmd, text=True)
    return json.loads(output)


def is_model_changed_in_pr(pr_number):
    # Get the list of files changed in the PR
    # Note: The diff has a maximum limit of 3000 files, but since most of the
    # PRs only change a few files, this should be sufficient. If there are more
    # than 300 files changed, we may need to use the GitHub API to get the
    # list of changed files instead.
    cmd = ["gh", "pr", "diff", str(pr_number), "--name-only"]
    output = subprocess.check_output(cmd, text=True)
    changed_files = output.splitlines()
    # Check if 'model.xml' is in the list of changed files
    return "model.xml" in changed_files


if __name__ == "__main__":
    run_tests_on_prs(pr_start=89)
    print("Test results saved to growth_match_summary.csv")
    print("You can plot the results using the provided plotting script.")
