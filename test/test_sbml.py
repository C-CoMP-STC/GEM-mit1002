import unittest

import cobra


class TestValidSBML(unittest.TestCase):
    def test_valid_sbml(self):
        # Validate the SBML file with COBRApy
        results = cobra.io.validate_sbml_model('model.xml')

        # Check that the SBML file is valid
        # By checking that there are no errors in the 'SBML_ERROR' key
        # of the second element of the results tuple. The first element
        # is the model itself.
        # A valid file should have 0 errors.
        self.assertEqual(0, len(results[1]['SBML_ERROR']))
