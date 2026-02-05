import cobra

# Load the E. coli core model
model = cobra.io.read_sbml_model("/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/ecoli/iJO1366.xml")

# Run normal FBA
fba_sol = model.optimize()
fba_fluxes = fba_sol.fluxes
with open("fba_fluxes.json", "w") as f:
    fba_fluxes.to_json(f)

# Save names of key reactions
biomass_rxn_name = "BIOMASS_Ec_iJO1366_core_53p95M"
glc_ex_name = "EX_glc__D_e"

# Define the CAFBA parameters
w_c = 0  # The proteime cost of the catabolic sector
w_i = 8.3e-4  # The specific proteome cost for each enzyme (in gh/mmol)
w_r = 0.169  # The cost of the ribosome sector (in h)
phi_max = 0.484  # The total fraction of the proteome available for growth-dependent sectors

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
for reaction in model.reactions:
    if reaction.id not in [biomass_rxn_name, glc_ex_name]:
        total_proteome_cost_expression += w_i * reaction.flux_expression

# Add the CAFBA constraint
cafba_constraint = model.problem.Constraint(
    total_proteome_cost_expression, lb=phi_max, ub=phi_max)
model.add_cons_vars(cafba_constraint)

# Optimize
cafba_sol = model.optimize()
cafba_fluxes = cafba_sol.fluxes
with open("cafba_fluxes.json", "w") as f:
    cafba_fluxes.to_json(f)
