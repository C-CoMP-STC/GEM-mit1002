import cobra
import numpy as np
import pandas as pd

# 1. Configuration Constants (from Paper/SI)
W_I = 8.3e-4  # Average metabolic weight
W_R = 0.17  # Ribosomal weight
PHI_Q = 0.45  # Housekeeping offset
PHI_R0 = 0.066  # Ribosomal offset
PHI_MAX = 1 - PHI_Q - PHI_R0  # Total allocatable budget (~0.484)

# Reaction IDs for my Amac model
BIOMASS_ID = "bio1_biomass"
GLC_ID = "EX_cpd00027_e0"
PFK_ID = "rxn00545_c0"  # EMP Pathway marker
EDD_ID = "rxn01477_c0"  # ED Pathway marker


# Helper function to run CAFBA with specified parameters
def run_cafba(
    model,
    w_c=0,
    w_i=8.3e-4,
    w_r=0.169,
    phi_max=0.484,
    biomass_id="BIOMASS_Ec_iJO1366_core_53p95M",
    glc_ex_id="EX_glc__D_e",
    ed_reactions=["rxn00604_c0", "rxn01476_c0", "rxn01477_c0", "rxn03884_c0"],
    ed_discount=3.5,
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
    - ed_reactions: A list of reaction IDs that belong to the ED pathway (default: ["G6PDH2r", "PGL", "EDD", "EDA"]).
    - ed_discount: The factor by which to discount the ED pathway enzymes' cost (default: 3.5).

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

        # Apply a discount to ED pathway enzymes
        if rxn.id in ed_reactions:
            w_i_effective = w_i / ed_discount
        else:
            w_i_effective = w_i

        # Add weight to both directions to handle |v|
        coefficients[rxn.forward_variable] = w_i_effective
        coefficients[rxn.reverse_variable] = w_i_effective

    # 4. Apply all coefficients at once
    cafba_cons.set_linear_coefficients(coefficients)

    # 5. Solve
    return model.optimize()


# --- Running the Carbon Limitation Simulation ---
# Load the model
model = cobra.io.read_sbml_model("../../model.xml")

# Set a minimal medium with glucose as the sole carbon source
minimal_media = {
    "EX_cpd00027_e0": 10,  # Glucose_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00007_e0": 20,  # O2_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00075_e0": 1000,  # Nitrite_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}
model.medium = minimal_media

# Generate a list of w_c values to test, with more resolution at the low end (0 to 1)
# Use the cubic mesh from the paper to get better resolution at the low end
Npoints = 100
w_vec = 0.01 * np.linspace(0, 1, Npoints) + 0.99 * np.power(
    np.linspace(0, 1, Npoints), 3
)

summary_data = []

print("Starting sweep...")
for wc in w_vec:
    sol = run_cafba(
        model, wc, biomass_id=BIOMASS_ID, glc_ex_id=GLC_ID, ed_reactions=[EDD_ID]
    )

    if sol.status == "optimal":
        # Calculate Sectors
        lambda_val = sol.objective_value
        phi_R = PHI_R0 + (W_R * lambda_val)
        phi_C = wc * abs(sol.fluxes[GLC_ID])

        # Calculate phi_E by summing weighted fluxes of all internal rxns
        internal_rxns = [
            r.id
            for r in model.reactions
            if r.id not in [BIOMASS_ID, GLC_ID]
            and not r.boundary
            and len(r.metabolites) >= 2
        ]
        phi_E = (sol.fluxes[internal_rxns].abs() * W_I).sum()

        # Save results
        summary_data.append(
            {
                "w_c": wc,
                "growth_rate": lambda_val,
                "phi_C": phi_C,
                "phi_E": phi_E,
                "phi_R": phi_R,
                "flux_PFK": sol.fluxes[PFK_ID],
                "flux_EDD": sol.fluxes[EDD_ID],
                "flux_glc": abs(sol.fluxes[GLC_ID]),
            }
        )

df_final = pd.DataFrame(summary_data)
df_final.to_csv("cafba_sweep_results.csv")
print("Done! Results saved to cafba_sweep_results.csv")
