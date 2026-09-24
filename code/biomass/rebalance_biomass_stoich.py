import os
import sys

import cobra
from gem_utilities.biomass import calculate_biomass_weight

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from tools.paths import MODEL_PATH  # noqa: E402

# Run the biomass weight calculation on the model
model = cobra.io.read_sbml_model(MODEL_PATH)
weight = calculate_biomass_weight(
    model,
    mets_to_ignore=["cpd11416_c0"],
)

# Rebalance the reaction stoichiometry
# NOTE: This also changes the coefficient of the biomass pseudometabolite,
# which it should NOT, that should be left as 1
# For now, I will just leave it like this, but you do need to go change just
# that coefficient back to 1 manually after running this script
# If we need to run this script multiple times, we should add code to reset
# the biomass pseudometabolite coefficient to 1 after rebalancing
if weight != 1.000:
    print(
        f"Biomass weight is {weight} g/mol, rebalancing biomass reaction stoichiometry"
    )
    biomass_reaction = model.reactions.get_by_id("bio1_biomass")
    biomass_reaction.add_metabolites(
        {
            met: coef * (1 / weight)
            for met, coef in biomass_reaction.metabolites.items()
        },
        combine=False,
    )

# Save the updated model
cobra.io.write_sbml_model(model, MODEL_PATH)
