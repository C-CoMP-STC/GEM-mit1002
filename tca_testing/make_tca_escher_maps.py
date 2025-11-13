import itertools

import cobra
import pandas as pd

import escher

# Load the models
amac_model = cobra.io.read_sbml_model("../model.xml")
ecoli_model = cobra.io.read_sbml_model(
    "/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/ecoli/iJO1366.xml"
)

# Save the paths to the maps
amac_map_path = "../escher/MIT1002_TCA_escher-map.json"
ecoli_map_path = (
    "/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/ecoli/iJO1366.TCA Only.json"
)

# Define media without carbon or nitrogen sources
amac_basal_media = {
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00007_e0": 20,  # O2_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}
ecoli_basal_media = {
    "EX_co2_e": 1000.0,
    "EX_cobalt2_e": 1000.0,
    "EX_h_e": 1000.0,
    "EX_h2o_e": 1000.0,
    "EX_k_e": 1000.0,
    "EX_cu2_e": 1000.0,
    "EX_mg2_e": 1000.0,
    "EX_mn2_e": 1000.0,
    "EX_mobd_e": 1000.0,
    "EX_na1_e": 1000.0,
    "EX_ca2_e": 1000.0,
    "EX_cbl1_e": 0.01,
    "EX_ni2_e": 1000.0,
    "EX_o2_e": 1000.0,
    "EX_cl_e": 1000.0,
    "EX_pi_e": 1000.0,
    "EX_zn2_e": 1000.0,
    "EX_sel_e": 1000.0,
    "EX_slnt_e": 1000.0,
    "EX_so4_e": 1000.0,
    "EX_tungs_e": 1000.0,
    "EX_fe2_e": 1000.0,
    "EX_fe3_e": 1000.0,
}

# Define additions to the minimal media
c_sources = {
    "Glucose": {"ModelSEED": "EX_cpd00027_e0", "BiGG": "EX_glc__D_e"},
    "Acetate": {"ModelSEED": "EX_cpd00029_e0", "BiGG": "EX_ac_e"},
}
n_sources = {
    "Ammonia": {"ModelSEED": "EX_cpd00013_e0", "BiGG": "EX_nh4_e"},
    "Nitrate": {"ModelSEED": "EX_cpd00209_e0", "BiGG": "EX_no3_e"},
}

# Make list of all reactions to maximize
modelseed_tca_rxns = [
    "rxn00256_c0",
    "rxn00974_c0",
    "rxn01388_c0",
    "rxn01387_c0",
    "rxn00199_c0",
    "rxn00336_c0",
    "rxn00330_c0",
    "rxn00441_c0",
    "rxn02376_c0",
    "rxn01872_c0",
    "rxn08094_c0",
    "rxn00285_c0",
    "rxn00288_c0",
    "rxn10126_c0",
    "rxn00799_c0",
    "rxn00248_c0",
]
bigg_tca_rxns = [
    "ACONTa",
    "ACONTb",
    "ICDHyr",
    "AKGDH",
    "SUCOAS",
    "FRD2",
    "SUCDi",
    "FUM",
    "MOX",
    "MDH2",
    "MDH",
    "MDH3",
    "CS",
    "CITL",
    "MALS",
    "ICL",
]

# Make interventions in the Amac model
# Lumped reactions for ICDH and AKGDH
# Make a copy of the model to work with
amac_w_lumped_rxns = amac_model.copy()

# Remove (not just knock out) current "unlumped" reactions
unlumped_reactions = [
    "rxn01387_c0",
    "rxn00199_c0",
    "rxn00441_c0",
    "rxn02376_c0",
    "rxn01872_c0",
]
amac_w_lumped_rxns.remove_reactions(unlumped_reactions)

# Define a lumped ICDH
lumped_icdh_reaction = cobra.Reaction("rxn00198_c0")
lumped_icdh_reaction.name = "NADP+-dependent isocitrate dehydrogenase"
lumped_icdh_reaction.lower_bound = 0  # irreversible
lumped_icdh_reaction.upper_bound = 1000  # arbitrary upper bound
lumped_icdh_reaction.add_metabolites(
    {
        amac_w_lumped_rxns.metabolites.cpd00260_c0: -1,  # isocitrate
        amac_w_lumped_rxns.metabolites.cpd00006_c0: -1,  # NADP
        amac_w_lumped_rxns.metabolites.cpd00024_c0: 1,  # alpha-ketoglutarate
        amac_w_lumped_rxns.metabolites.cpd00011_c0: 1,  # CO2
        amac_w_lumped_rxns.metabolites.cpd00005_c0: 1,  # NADPH
    }
)

# Define a lumped AKGDH
lumped_akgdh_rxn = cobra.Reaction("rxn08094_c0")
lumped_akgdh_rxn.name = "2-Oxoglutarate dehydrogenase complex"
lumped_akgdh_rxn.lower_bound = 0
lumped_akgdh_rxn.upper_bound = 1000
lumped_akgdh_rxn.add_metabolites(
    {
        amac_w_lumped_rxns.metabolites.cpd00024_c0: -1.0,  # 2-Oxoglutarate
        amac_w_lumped_rxns.metabolites.cpd00010_c0: -1.0,  # CoA
        amac_w_lumped_rxns.metabolites.cpd00003_c0: -1.0,  # NAD+
        amac_w_lumped_rxns.metabolites.cpd00078_c0: 1.0,  # Succinyl-CoA
        amac_w_lumped_rxns.metabolites.cpd00011_c0: 1.0,  # CO2
        amac_w_lumped_rxns.metabolites.cpd00004_c0: 1.0,  # NADH
    }
)

# Add lumped reactions
amac_w_lumped_rxns.add_reactions([lumped_icdh_reaction, lumped_akgdh_rxn])

# TODO: Strict ATP Production
# TODO: List reactions involving ATP to make irreversible (force in the ATP-consuming direction)
atp_consuming_reactions = []

# Optimize for Biomass (Check fluxes when growing)
# For Amac
# Make a list of all versions of the amac model to test
amac_models_to_test = {
    "Original": amac_model,
    "With_Lumped_Reactions": amac_w_lumped_rxns,
}
# Make a list to store results for each model
amac_growth_dfs = {}
c_names = list(c_sources.keys())
n_names = list(n_sources.keys())

for model_name, model in amac_models_to_test.items():
    amac_growth_results = []
    print(f"Running AMAC tests for model: {model_name}")
    for c_k in range(1, len(c_names) + 1):
        for c_subset in itertools.combinations(c_names, c_k):
            for n_k in range(1, len(n_names) + 1):
                for n_subset in itertools.combinations(n_names, n_k):
                    c_label = "+".join(c_subset)
                    n_label = "+".join(n_subset)
                    print(f"  Testing C source(s): {c_label} | N source(s): {n_label}")
                    # Make a medium (copy so basal isn't mutated)
                    amac_medium = amac_basal_media.copy()
                    # Add carbon sources (use same per-source bounds as before)
                    for cname in c_subset:
                        cid = c_sources[cname]["ModelSEED"]
                        amac_medium[cid] = 10 if cname == "Glucose" else 30
                    # Add nitrogen sources (unlimited)
                    for nname in n_subset:
                        nid = n_sources[nname]["ModelSEED"]
                        amac_medium[nid] = 1000
                    # Work on a fresh copy of the model to avoid mutating originals
                    m = model.copy()
                    m.medium = amac_medium
                    sol = m.optimize()
                    # Collect TCA fluxes
                    rxn_fluxes = {
                        rxn: sol.fluxes[rxn]
                        for rxn in modelseed_tca_rxns
                        if rxn in sol.fluxes
                    }
                    rxn_fluxes["bio1_biomass"] = sol.objective_value
                    amac_growth_results.append(
                        {
                            "Model": model_name,
                            "C_source": c_label,
                            "N_source": n_label,
                            "fluxes": rxn_fluxes,
                        }
                    )
    # Convert to DataFrame for this model
    amac_growth_dfs[model_name] = pd.DataFrame(amac_growth_results)

# Make a save all escher maps as HTML files
for model_name, model in amac_models_to_test.items():
    growth_df = amac_growth_dfs[model_name]
    for c_k in range(1, len(c_names) + 1):
        for c_subset in itertools.combinations(c_names, c_k):
            for n_k in range(1, len(n_names) + 1):
                for n_subset in itertools.combinations(n_names, n_k):
                    c_label = "+".join(c_subset)
                    n_label = "+".join(n_subset)
                    flux_data = growth_df[
                        (growth_df["C_source"] == c_label)
                        & (growth_df["N_source"] == n_label)
                    ]["fluxes"].item()
                    escher.Builder(
                        model=model, map_json=amac_map_path, reaction_data=flux_data
                    ).save_html(
                        f"escher_plots/html/amac_{model_name}_{c_label}+{n_label}.html"
                    )

# For e. coli
# Make a list to store results
ecoli_growth_results = []

# Iterate over all non-empty subsets of carbon and nitrogen sources
for c_k in range(1, len(c_names) + 1):
    for c_subset in itertools.combinations(c_names, c_k):
        for n_k in range(1, len(n_names) + 1):
            for n_subset in itertools.combinations(n_names, n_k):
                c_label = "+".join(c_subset)
                n_label = "+".join(n_subset)
                # Debugging
                print(
                    f"Testing E. coli with C source: {c_label} and N source: {n_label}"
                )
                # Make a medium
                ecoli_medium = ecoli_basal_media.copy()
                # Add carbon sources (use same per-source bounds as before)
                for cname in c_subset:
                    cid = c_sources[cname]["BiGG"]
                    if cname == "Glucose":
                        ecoli_medium[cid] = 10
                    else:
                        ecoli_medium[cid] = 30
                # Add nitrogen sources (unlimited)
                for nname in n_subset:
                    nid = n_sources[nname]["BiGG"]
                    ecoli_medium[nid] = 1000
                # Set the medium
                ecoli_model.medium = ecoli_medium
                # Optimize for Biomass (Check Fluxes When Growing)
                ecoli_solution = ecoli_model.optimize()
                # Get the fluxes through the TCA cycle reactions (skip reactions that aren't in the model)
                rxn_fluxes = {
                    rxn: ecoli_solution.fluxes[rxn]
                    for rxn in bigg_tca_rxns
                    if rxn in ecoli_solution.fluxes
                }
                # Add the biomass flux to the dict
                rxn_fluxes["BIOMASS_Ec_iJO1366_core_53p95M"] = (
                    ecoli_solution.objective_value
                )
                # Store the results
                ecoli_growth_results.append(
                    {
                        "C_source": c_label,
                        "N_source": n_label,
                        "fluxes": rxn_fluxes,
                    }
                )
# Convert to DataFrame for easier viewing
ecoli_growth_df = pd.DataFrame(ecoli_growth_results)

# Make a save all escher maps as HTML files
for c_k in range(1, len(c_names) + 1):
    for c_subset in itertools.combinations(c_names, c_k):
        for n_k in range(1, len(n_names) + 1):
            for n_subset in itertools.combinations(n_names, n_k):
                c_label = "+".join(c_subset)
                n_label = "+".join(n_subset)
                flux_data = ecoli_growth_df[
                    (ecoli_growth_df["C_source"] == c_label)
                    & (ecoli_growth_df["N_source"] == n_label)
                ]["fluxes"].item()
                escher.Builder(
                    model=ecoli_model, map_json=ecoli_map_path, reaction_data=flux_data
                ).save_html(f"escher_plots/html/ecoli_{c_label}+{n_label}.html")

# Optimize for TCA cycle reactions (check for blockages)
# Make a list to store results
amac_blockage_results = []
# Get every combination of C and N source
for c_name, c_id_dict in c_sources.items():
    for n_name, n_id_dict in n_sources.items():
        # Debugging
        print(
            f'Testing AMAC with C source: {c_name} ({c_id_dict["ModelSEED"]}) and N source: {n_name} ({n_id_dict["ModelSEED"]})'
        )
        # Make a medium
        amac_medium = amac_basal_media.copy()
        # Add a limiting amount of carbon source
        # But use different bounds for glucose and acetate so the number of carbon atoms is the same
        if c_name == "Glucose":
            amac_medium[c_id_dict["ModelSEED"]] = 10
        else:
            amac_medium[c_id_dict["ModelSEED"]] = 30
        amac_medium[n_id_dict["ModelSEED"]] = 1000
        # Set the medium
        amac_model.medium = amac_medium
        # Make a dictionary to store individual reaction fluxes
        rxn_fluxes = {}
        # Optimize for each of the TCA cycle reactions individually
        for rxn_id in modelseed_tca_rxns:
            # Skip if the reaction is not in the model
            if rxn_id not in [r.id for r in amac_model.reactions]:
                continue
            # Set the objective to the current reaction
            amac_model.objective = rxn_id
            # Debugging
            print(f"  Optimizing for reaction: {rxn_id}")
            # Optimize
            amac_solution = amac_model.optimize()
            # Store the flux
            rxn_fluxes[rxn_id] = amac_solution.objective_value
        # Store the results
        amac_blockage_results.append(
            {
                "C_source": c_name,
                "N_source": n_name,
                "fluxes": rxn_fluxes,
            }
        )
# Convert to DataFrame for easier viewing
amac_blockage_df = pd.DataFrame(amac_blockage_results)

# Make a save all escher maps as HTML files
for c_label in c_names:
    for n_label in n_names:
        flux_data = amac_blockage_df[
            (amac_blockage_df["C_source"] == c_label)
            & (amac_blockage_df["N_source"] == n_label)
        ]["fluxes"].item()
        escher.Builder(
            model=amac_model, map_json=amac_map_path, reaction_data=flux_data
        ).save_html(
            f"escher_plots/html/amac_Original_{c_label}+{n_label}_blocked_reactions.html"
        )

# For E. coli
# Make a list to store results
ecoli_blockage_results = []
# Get every combination of C and N source
for c_name, c_id_dict in c_sources.items():
    for n_name, n_id_dict in n_sources.items():
        # Debugging
        print(
            f'Testing E. coli with C source: {c_name} ({c_id_dict["BiGG"]}) and N source: {n_name} ({n_id_dict["BiGG"]})'
        )
        # Make a medium
        ecoli_medium = ecoli_basal_media.copy()
        # Add a limiting amount of carbon source
        # But use different bounds for glucose and acetate so the number of carbon atoms is the same
        if c_name == "Glucose":
            ecoli_medium[c_id_dict["BiGG"]] = 10
        else:
            ecoli_medium[c_id_dict["BiGG"]] = 30
        # Add an unlimited amount of nitrogen source
        ecoli_medium[n_id_dict["BiGG"]] = 1000
        # Set the medium
        ecoli_model.medium = ecoli_medium
        # Make a dictionary to store individual reaction fluxes
        rxn_fluxes = {}
        # Optimize for each of the TCA cycle reactions individually
        for rxn_id in bigg_tca_rxns:
            # Skip if the reaction is not in the model
            if rxn_id not in [r.id for r in ecoli_model.reactions]:
                continue
            # Set the objective to the current reaction
            ecoli_model.objective = rxn_id
            # Debugging
            print(f"  Optimizing for reaction: {rxn_id}")
            # Optimize
            ecoli_solution = ecoli_model.optimize()
            # Store the flux
            rxn_fluxes[rxn_id] = ecoli_solution.objective_value
        # Store the results
        ecoli_blockage_results.append(
            {
                "C_source": c_name,
                "N_source": n_name,
                "fluxes": rxn_fluxes,
            }
        )
# Convert to DataFrame for easier viewing
ecoli_blockage_df = pd.DataFrame(ecoli_blockage_results)

# Make a save all escher maps as HTML files
for c_label in c_names:
    for n_label in n_names:
        flux_data = ecoli_blockage_df[
            (ecoli_blockage_df["C_source"] == c_label)
            & (ecoli_blockage_df["N_source"] == n_label)
        ]["fluxes"].item()
        escher.Builder(
            model=ecoli_model, map_json=ecoli_map_path, reaction_data=flux_data
        ).save_html(
            f"escher_plots/html/ecoli_{c_label}+{n_label}_blocked_reactions.html"
        )
