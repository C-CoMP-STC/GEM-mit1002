import cobra
import numpy as np
import pandas as pd


# Helper function to run CAFBA with specified parameters
def run_cafba(
    model,
    w_c=0,
    w_i=8.3e-4,
    w_r=0.169,
    phi_max=0.484,
    biomass_id="BIOMASS_Ec_iJO1366_core_53p95M",
    glc_ex_id="EX_glc__D_e",
):
    """
    Run CAFBA on the given model with the specified parameters.

    Parameters:
    - model: A COBRApy model object.
    - w_c: The proteome cost of the catabolic sector (default: 0).
    - w_i: The specific proteome cost for each enzyme (in gh/mmol, default: 8.3e-4).
    - w_r: The cost of the ribosome sector (in h, default: 0.169).
    - phi_max: The total fraction of the proteome available for growth-dependent sectors (default: 0.484).
    - biomass_id: The ID of the biomass reaction in the model (default: "BIOMASS_Ec_iJO1366_core_53p95M").
    - glc_ex_id: The ID of the glucose exchange reaction in the model (default: "EX_glc__D_e").

    Returns:
    A COBRA solution containing the results from the CAFBA solution.
    """

    # 1. Clear existing CAFBA constraints if they exist
    if "CAFBA_constraint" in model.constraints:
        model.remove_cons_vars(model.constraints.CAFBA_constraint)

    # 2. Define the constraint (initially empty)
    # We use an inequality: Sum of costs <= phi_max
    cafba_cons = model.problem.Constraint(0, lb=0, ub=phi_max, name="CAFBA_constraint")
    model.add_cons_vars(cafba_cons)

    # 3. Build the coefficients dictionary (The FAST way)
    # This maps the internal solver variables to their weights
    coefficients = {}

    # Sector R: Ribosomes (Weight the biomass flux)
    biomass_rxn = model.reactions.get_by_id(biomass_id)
    coefficients[biomass_rxn.forward_variable] = w_r

    # Sector C: Carbon uptake (Weight the absolute value of glucose uptake)
    # Since uptake is negative, we weight the reverse_variable (the 'in' direction)
    glc_rxn = model.reactions.get_by_id(glc_ex_id)
    coefficients[glc_rxn.reverse_variable] = w_c

    # Sector E: Internal Enzymes
    for rxn in model.reactions:
        # Avoid double-counting biomass and the specific carbon source
        if rxn.id == biomass_id or rxn.id == glc_ex_id:
            continue

        # Ignore exchange/sink/demand reactions (they don't have enzymes)
        if len(rxn.metabolites) < 2 or rxn.boundary:
            continue

        # Add weight to both directions to handle |v|
        coefficients[rxn.forward_variable] = w_i
        coefficients[rxn.reverse_variable] = w_i

    # 4. Apply all coefficients at once
    cafba_cons.set_linear_coefficients(coefficients)

    # 5. Solve
    return model.optimize()


# --- Running the Carbon Limitation Simulation ---
# Load the E. coli model
model = cobra.io.load_json_model(
    "/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/ecoli/iJO1366.json"
)
# Generate a list of w_c values to test, with more resolution at the low end (0 to 1)
# Use the cubic mesh from the paper to get better resolution at the low end
Npoints = 100
w_vec = 0.01 * np.linspace(0, 1, Npoints) + 0.99 * np.power(
    np.linspace(0, 1, Npoints), 3
)
# Make a dictionary to store the results for different wc values
results = {}
# Loop through the values of wc and run CAFBA for each value, storing the results in the dictionary
for wc in w_vec:
    print(f"Running CAFBA with wc={wc:.2f}...")
    cafba_solution = run_cafba(model, w_c=wc)
    results[wc] = cafba_solution.fluxes
# Convert the results dictionary to a DataFrame for easier analysis and visualization
results_df = pd.DataFrame(results)
# Save the results to a CSV file
results_df.to_csv("cafba_results.csv")
