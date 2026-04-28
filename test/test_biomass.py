import unittest

import cobra
from gem_utilities.biomass import calculate_biomass_weight


class TestBiomass(unittest.TestCase):
    def test_biomass_weight(self):
        """Test that the model's biomass metabolite is defined to be 1 g"""
        # Load the model
        model = cobra.io.read_sbml_model("model.xml")

        # Calculate the biomass weight
        biomass_weight = calculate_biomass_weight(
            model=model, mets_to_ignore=["cpd11416_c0"]
        )

        # Assert that the biomass weight is approximately 1 g
        self.assertAlmostEqual(
            biomass_weight,
            1.0,
            places=6,
            msg=f"Biomass weight is {biomass_weight}, not 1 g",
        )
