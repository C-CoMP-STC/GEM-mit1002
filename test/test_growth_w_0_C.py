import unittest

import cobra

class TestGrowthWithNoCarbon(unittest.TestCase):
    def test_growth_w_0_C(self):
        # Load the model with cobrapy
        model = cobra.io.read_sbml_model('../model.xml')

        # Set the media so that there are no carbon sources
        medium = {'EX_co2_e': 1000.0,
                  'EX_glc__D_e': 0.0,
                  'EX_h_e': 1000.0,
                  'EX_h2o_e': 1000.0,
                  'EX_nh4_e': 1000.0,
                  'EX_o2_e': 1000.0,
                  'EX_pi_e': 1000.0}

        model.medium = medium

        # Run the model
        sol = model.optimize()

        # Check that no biomass is produced
        self.assertEqual(sol.objective_value, 0)

if __name__ == '__main__':
    unittest.main()
