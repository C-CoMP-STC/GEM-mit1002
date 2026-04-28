import unittest

import cobra


class TestValidSBML(unittest.TestCase):
    def test_valid_sbml(self):
        # Validate the SBML file with COBRApy
        results = cobra.io.validate_sbml_model("model.xml")

        # Check that the SBML file is valid
        # By checking that there are no errors in the 'SBML_ERROR' key
        # of the second element of the results tuple. The first element
        # is the model itself.
        # A valid file should have 0 errors.
        errors = results[1]["SBML_ERROR"]
        self.assertEqual(0, len(errors), msg=f"SBML validation errors: {errors}")

    def test_isolated_genes_and_mets(self):
        # Load the model
        model = cobra.io.read_sbml_model("model.xml")

        # Check for isolated genes
        isolated_genes = [g.id for g in model.genes if len(g.reactions) == 0]
        self.assertEqual(
            len(isolated_genes),
            0,
            msg=f"{len(isolated_genes)} isolated gene(s) found: {isolated_genes}",
        )

        # Check for isolated metabolites
        isolated_mets = [m.id for m in model.metabolites if len(m.reactions) == 0]
        self.assertEqual(
            len(isolated_mets),
            0,
            msg=f"{len(isolated_mets)} isolated metabolite(s) found: {isolated_mets}",
        )

    def test_mass_balance(self):
        # Load the model
        model = cobra.io.read_sbml_model("model.xml")

        # Check for mass balance in each reaction
        results = cobra.manipulation.check_mass_balance(model)

        # Check that there are no reactions returned by the mass balance check
        self.assertEqual(
            len(results),
            0,
            msg=f"{len(results)} reactions are not mass balanced: {results}",
        )
