import requests
from cobra import Reaction, Metabolite


def fetch_kegg_reaction(kegg_id):
    # KEGG API URL for fetching reaction details
    url = f"https://rest.kegg.jp/get/rn:{kegg_id}"
    response = requests.get(url)
    if response.status_code != 200:
        print("Failed to fetch data from KEGG API")
        return None

    # Parse the response text to find the reaction equation
    lines = response.text.split("\n")
    equation = None
    for line in lines:
        if line.startswith("EQUATION"):
            equation = line.split(maxsplit=1)[1]
            break

    if not equation:
        print("Equation not found in KEGG response")
        return None

    # Create a COBRApy reaction object
    reaction = Reaction(kegg_id)
    reaction.name = kegg_id
    reaction.subsystem = "KEGG Reaction"
    reaction.lower_bound = (
        0  # This might change depending on the reaction reversibility
    )
    reaction.upper_bound = 1000.0

    # Parse the equation for reactants and products
    left_side, right_side = equation.split(" <=> ")
    for metabolite_str in left_side.split(" + "):
        coeff, met_id = metabolite_str.split()
        coeff = -float(coeff)  # Reactants have negative stoichiometry
        metabolite = Metabolite(
            met_id.strip(), formula="", name=met_id.strip(), compartment="c"
        )
        reaction.add_metabolites({metabolite: coeff})

    for metabolite_str in right_side.split(" + "):
        coeff, met_id = metabolite_str.split()
        coeff = float(coeff)
        metabolite = Metabolite(
            met_id.strip(), formula="", name=met_id.strip(), compartment="c"
        )
        reaction.add_metabolites({metabolite: coeff})

    return reaction


# Example usage
kegg_reaction_id = "R00001"
reaction = fetch_kegg_reaction(kegg_reaction_id)
if reaction:
    print(f"Reaction ID: {reaction.id}")
    print(f"Reaction: {reaction.build_reaction_string(use_metabolite_names=False)}")
