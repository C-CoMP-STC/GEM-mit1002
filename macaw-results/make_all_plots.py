import os

import matplotlib.pyplot as plt
import pandas as pd

# Set the output directory (where the results.pkl file will be saved)
OUT_DIR = os.path.dirname(os.path.realpath(__file__))

# Set a folder for the plots
OUTPUT_FOLDER = os.path.join(OUT_DIR, "plots")


def main():
    # Check if the output folder exists, if not, create it
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    # Load in the results csv
    results = pd.read_csv(os.path.join(OUT_DIR, "helen-alteromonas_test-results.csv"))

    # Make a barchart with the proportion of reactions that fail each test (has any value other than "ok")

    # For the test columns (columns 'duplicate_test_exact',
    # 'duplicate_test_directions', 'duplicate_test_coefficients',
    # 'duplicate_test_redox', 'diphosphate_test', 'loop_test', 'dilution_test',
    # 'dead_end_test', 'pathway') calculate the proportion of reactions that
    # fail the test (have any value other than "ok" in the column) and divide
    # by the number of reactions in the column
    test_results = results[["duplicate_test_exact", "duplicate_test_directions", "duplicate_test_coefficients",
                            "duplicate_test_redox", "diphosphate_test", "loop_test", "dilution_test", "dead_end_test"]]
    test_results = test_results.applymap(lambda x: 1 if x != "ok" else 0)
    test_results = test_results.sum()/len(results)

    # Plot
    test_results.plot(kind="bar")
    plt.ylabel("Proportion of reactions failing test")
    plt.title("MACAW Test Results")
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_FOLDER, "test_results.png"))


if __name__ == "__main__":
    main()
