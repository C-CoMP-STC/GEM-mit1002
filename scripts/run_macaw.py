# This script is run in the custom CI workflow, and the output is saved

import cobra
from macaw.main import run_all_tests

# import py4cytoscape as p4c


# Run the MACAW pipeline
def run_macaw(model):
    # TODO: Figure out how to run the function without having to wrap it
    # Run all tests
    (test_results, edge_list) = run_all_tests(model)

    return test_results, edge_list


if __name__ == "__main__":
    # Load the model
    model = cobra.io.read_sbml_model("model.xml")

    # Run the MACAW pipeline
    (test_results, edge_list) = run_all_tests(model)

    # Save the results
    test_results.to_csv("macaw_results.csv")
    edge_list.to_csv("macaw_edge_list.csv")

    # TODO: Visualize the network results
