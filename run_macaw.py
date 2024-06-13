# This script is run in the custom CI workflow, and the output is saved

import cobra
from macaw import macaw_main

# Load the model
model = cobra.io.read_sbml_model("model.xml")

# Run the MACAW pipeline
(test_results, edge_list) = macaw_main.run_all_tests(model)

# Save the results
test_results.to_csv("macaw_results.csv")
edge_list.to_csv("macaw_edge_list.csv")
