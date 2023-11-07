import unittest

import cobra
import memote.support.consistency as consistency


class TestCycles(unittest.TestCase):
    # Check that there are no ATP generating cycles
    def test_atp_generating_cycles(self):
        # Load the model
        model = cobra.io.read_sbml_model("model.xml")

        # Set the metabolite to ATP
        met = "MNXM3"  # ATP

        # Run just the function for finding ATP generating cycles from MEMOTE
        # Don't generate the whole report
        rxns = consistency.detect_energy_generating_cycles(model, met)

        self.assertEqual(
            len(rxns),
            0,
            "Expect there to be 0 reactions involved in ATP generating cycles, but there are "
            + str(len(rxns)),
        )


if __name__ == "__main__":
    unittest.main()
