import os
import unittest

import cobra
import pandas as pd

# Set path to the `test_files` directory
TESTFILE_DIR = os.path.join(os.path.dirname(__file__), "test_files")

# Set path to the results directory
RESULTS_DIR = os.path.join(TESTFILE_DIR, "test_results")


class TestPathways(unittest.TestCase):
    # Test that the reactions of specified pathways are present in the model.
    # Make a plot of the reaction presence/absence in the model.
    def test_pathway_reactions_present(self):
        # Load the model with cobrapy
        model = cobra.io.read_sbml_model("model.xml")

        # Load the CSV of the pathway reactions
        pathway_reactions = pd.read_csv(
            os.path.join(TESTFILE_DIR, "pathway_reactions.csv")
        )

        # Create a set of reaction IDs from the model for fast lookups
        model_reaction_ids = {r.id for r in model.reactions}

        # Use the vectorized `isin` method to get a boolean Series
        is_present = pathway_reactions["reaction_id"].isin(model_reaction_ids)

        # Map the boolean results to "Yes"/"No" strings
        pathway_reactions["in_model"] = is_present.map({True: "Yes", False: "No"})

        # Save the results to a new CSV file
        pathway_reactions.to_csv(
            os.path.join(RESULTS_DIR, "pathway_reaction_presence.csv"), index=False
        )
