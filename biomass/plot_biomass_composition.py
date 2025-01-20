import os

import cobra
import matplotlib.pyplot as plt

FILEPATH = os.path.dirname(__file__)

# Define the file path, and make it if it doesn't already exist
outpath = os.path.join(FILEPATH, "plots")
if not os.path.exists(outpath):
    os.makedirs(outpath)


# Write function to plot the biomass compostion as a pie chart
def plot_biomass_composition(model, biomass_rxn_id='bio1_biomass', outpath=outpath, filename="biomass_composition.png"):
    # Get the biomass reaction
    biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)

    # Get the biomass composition
    biomass_composition = biomass_rxn.metabolites

    # Compute mass contributions
    # FIXME: Why don't these add to 1?
    input_contributions = {}
    output_contributions = {}
    for metabolite, coeff in biomass_composition.items():
        molecular_weight = metabolite.formula_weight  # Get molecular weight
        # Get the mass contribution of the metabolite
        met_mass = abs(coeff) * molecular_weight
        if coeff < 0:
            input_contributions[metabolite.name] = met_mass
        else:
            output_contributions[metabolite.name] = met_mass

    # Normalize to percentages
    # TO DO: Make this a spearate function, so I can call it for the two different dictionaries
    total_mass = sum(input_contributions.values())
    # I checked that this does indeed add up to 100
    mass_percentages = {met: (mass / total_mass) * 100 for met, mass in input_contributions.items()}

    # Sort and prepare data for plotting
    labels, sizes = zip(*sorted(mass_percentages.items(), key=lambda x: x[1], reverse=True))

    # Create pie chart
    fig, ax = plt.subplots(figsize=(10, 6))
    wedges, texts = ax.pie(sizes, startangle=140)

    # Add a legend with metabolite names and corresponding colors
    ax.legend(wedges, labels, title="Biomass Components", loc="center left", bbox_to_anchor=(1, 0.5))

    plt.title("Biomass Composition (Mass %)")
    plt.savefig(os.path.join(outpath, filename))


# Test: Plot the biomass composition of the E. coli core model
model = cobra.io.load_model("textbook")
plot_biomass_composition(model, "Biomass_Ecoli_core", filename="biomass_composition_ecoli_core.png")

# Run it on the KBase model
model = cobra.io.read_sbml_model("/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/mit1002-model/2025-01-08_Scott_draft-model-from-KBase.xml")
plot_biomass_composition(model, filename="biomass_composition_kbase.png")
