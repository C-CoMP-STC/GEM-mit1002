import os
import sys

import cobra
from gem_utilities.biomass import save_biomass_composition_work_table

# Make `tools` importable; everything else comes from tools.paths.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tools.paths import MODEL_PATH  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "results")

# Load the model
model = cobra.io.read_sbml_model(MODEL_PATH)

# Save the biomass composition table
save_biomass_composition_work_table(
    model=model, mets_to_ignore=["cpd11416_c0"], out_dir=RESULTS_DIR
)
