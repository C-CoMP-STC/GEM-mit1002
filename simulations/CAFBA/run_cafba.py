import cobra


def run_cafba(
    model,
    w_c=0,
    w_i=8.3e-4,
    w_r=0.169,
    phi_max=0.484,
    biomass_rxn_name="BIOMASS_Ec_iJO1366_core_53p95M",
    glc_ex_name="EX_glc__D_e",
):
    """
    Run CAFBA on the given model with the specified parameters.

    Parameters:
    - model: A COBRApy model object.
    - w_c: The proteome cost of the catabolic sector (default: 0).
    - w_i: The specific proteome cost for each enzyme (in gh/mmol, default: 8.3e-4).
    - w_r: The cost of the ribosome sector (in h, default: 0.169).
    - phi_max: The total fraction of the proteome available for growth-dependent sectors (default: 0.484).

    Returns:
    - cafba_fluxes: A pandas Series containing the fluxes from the CAFBA solution.
    """

    # Create a variable to hold all of the proteome constraint
    total_proteome_cost_expression = 0

    # Add the ribosomal cost
    biomass_rxn = model.reactions.get_by_id(biomass_rxn_name)
    total_proteome_cost_expression += w_r * biomass_rxn.flux_expression

    # Add the catabolic cost
    carbon_uptake_rxn = model.reactions.get_by_id(glc_ex_name)
    total_proteome_cost_expression += w_c * carbon_uptake_rxn.flux_expression

    # Add enzyme costs for all other reactions
    # TODO: Should I exclude certain reactions here (e.g., exchange reactions, biomass reaction, etc.)?
    # TODO: Should I use subsystem or is_transport to exclude certain reactions?
    for reaction in model.reactions:
        if reaction.id not in [biomass_rxn_name, glc_ex_name]:
            total_proteome_cost_expression += w_i * reaction.flux_expression

    # Add the CAFBA constraint
    cafba_constraint = model.problem.Constraint(
        total_proteome_cost_expression, lb=phi_max, ub=phi_max
    )
    model.add_cons_vars(cafba_constraint)

    # Optimize
    cafba_sol = model.optimize()

    # Return the solution
    return cafba_sol


# Load the E. coli model
model = cobra.io.load_json_model(
    "/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/ecoli/iJO1366.json"
)
