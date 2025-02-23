#!~/opt/miniconda3/envs/med4-hot1a3/bin/python

import os
import re
from collections import defaultdict

import pandas as pd

FILE_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    # Replace 'biomass_data.csv' with the path to your CSV file.
    data_file = os.path.join(FILE_DIR, "gram_negative_biomass_components.tsv")
    check_biomass_balance(data_file)


def parse_formula(formula):
    """
    Parses a chemical formula into a dictionary of element counts.
    For example, 'C6H12O6' becomes {'C': 6, 'H': 12, 'O': 6}.
    """
    pattern = r"([A-Z][a-z]?)(\d*)"
    counts = defaultdict(int)
    for element, count in re.findall(pattern, formula):
        counts[element] += int(count) if count else 1
    return dict(counts)


def check_biomass_balance(data_file):
    """
    Loads the biomass data from a CSV file and calculates the net elemental balance
    and net charge of the biomass reaction.
    The CSV file must have columns: component, formula, mass, charge, coefficient.
    """
    df = pd.read_csv(data_file, sep="\t", index_col=0)

    net_elements = defaultdict(float)
    net_charge = 0.0

    # Process each biomass component.
    for idx, row in df.iterrows():
        formula = row["formula"]
        coeff = row["coefficient"]
        charge = row["charge"]

        # If the formula is missing, skip this component.
        if pd.isnull(formula):
            print(f"Skipping elements fot component '{idx}' because it is missing a formula.")
        else:
            # Parse the formula into element counts.
            elem_counts = parse_formula(formula)

            # Multiply each element count by the stoichiometric coefficient.
            for element, count in elem_counts.items():
                net_elements[element] += coeff * count

        # If the charge is missing, skip this component.
        if pd.isnull(charge):
            print(f"Skipping charge for component '{idx}' because it is missing a charge.")
        else:
            # Sum the net charge contribution.
            net_charge += coeff * charge

    # Print the net elemental balances.
    print(
        "Net elemental balance (should be near zero for a balanced biomass reaction):"
    )
    for element, net_count in net_elements.items():
        print(f"  {element}: {net_count:.4f}")

    # Print the net charge.
    print(f"\nNet charge: {net_charge:.4f}")


if __name__ == "__main__":
    main()
