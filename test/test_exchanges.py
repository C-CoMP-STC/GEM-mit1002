import unittest

import cobra


class TestExchanges(unittest.TestCase):
    # Check that there are no dead-end transporters (i.e. external metabolites without an exchange reaction)
    def test_dead_end_transporters(self):
        # Load the model
        model = cobra.io.read_sbml_model("model.xml")

        # Find all external metabolites
        external_metabolites = [
            met for met in model.metabolites if met.compartment == "e0"
        ]

        # Find all exchange reactions
        exchange_reactions = [
            rxn for rxn in model.reactions if rxn.id.startswith("EX_")
        ]

        # Create a set of metabolites involved in exchange reactions
        exchange_metabolites = set()
        for rxn in exchange_reactions:
            for met in rxn.metabolites:
                exchange_metabolites.add(met.id)

        # Check for external metabolites without an exchange reaction
        dead_end_metabolites = [
            met.id for met in external_metabolites if met.id not in exchange_metabolites
        ]

        # Assert that there are no dead-end metabolites
        self.assertEqual(
            len(dead_end_metabolites),
            0,
            msg=f"Dead-end transporters found for metabolites: {dead_end_metabolites}",
        )


if __name__ == "__main__":
    unittest.main()
