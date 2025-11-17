#!/Users/helenscott/opt/miniconda3/envs/escher/bin/python3
"""
This script generates escher maps with fluxes while growing, and while
optimizing for the TCA cycle reactions individually, for the mit1002 model,
the model with some interventions to test, and the E. coli iJO1366 model.
"""
import itertools
import os

import cobra
import pandas as pd
from playwright.sync_api import sync_playwright

import escher

# Define constants and paths at the module level
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(FILE_DIR)
ESCHER_PLOT_DIR = os.path.join(FILE_DIR, "escher_plots")
HTML_DIR = os.path.join(ESCHER_PLOT_DIR, "html")
SVG_DIR = os.path.join(ESCHER_PLOT_DIR, "svg")
E_COLI_MODEL_PATH = "/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/ecoli"


def main():
    """
    Main function to run all simulations and generate Escher maps.
    """
    ####################################################################
    # SET UP
    ####################################################################
    # Create directory for escher plots if it doesn't exist
    os.makedirs(ESCHER_PLOT_DIR, exist_ok=True)
    os.makedirs(HTML_DIR, exist_ok=True)
    os.makedirs(SVG_DIR, exist_ok=True)

    # Load models and define map paths
    amac_model = cobra.io.read_sbml_model(os.path.join(PROJECT_ROOT, "model.xml"))
    ecoli_model = cobra.io.read_sbml_model(
        os.path.join(E_COLI_MODEL_PATH, "iJO1366.xml")
    )
    amac_map_path = os.path.join(PROJECT_ROOT, "escher", "MIT1002_TCA_escher-map.json")
    ecoli_map_path = os.path.join(FILE_DIR, "iJO1366_tca-only-escher-map.json")

    # Define media and reactions
    # See helper function at the bottom of the script
    (
        amac_basal_media,
        ecoli_basal_media,
        c_sources,
        n_sources,
        modelseed_tca_rxns,
        bigg_tca_rxns,
    ) = define_media_and_reactions()

    ####################################################################
    # MAKE INTERVENTIONS IN THE AMAC MODEL
    ####################################################################
    # --- Make a version of the model with lumped ICDH and AKGDH reactions ---
    # Make a copy of the original model to modify
    amac_w_lumped_rxns = amac_model.copy()
    # Remove (not just knockout) the original ICDH and AKGDH reactions
    amac_w_lumped_rxns.remove_reactions(
        ["rxn01387_c0", "rxn00199_c0", "rxn00441_c0", "rxn02376_c0", "rxn01872_c0"]
    )

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
    amac_w_lumped_rxns.add_reactions([lumped_icdh_reaction, lumped_akgdh_rxn])

    # --- Make a version of the model with strict ATP production ---
    # List reactions involving ATP to make irreversible (consumption only)
    atp_consuming_reactions = [
        "rxn00077_c0",
        "rxn00104_c0",
        "rxn00239_c0",
        "rxn00364_c0",
        "rxn00379_c0",
        "rxn01219_c0",
        "rxn01509_c0",
        "rxn01517_c0",
        "rxn02314_c0",
        "rxn08762_c0",
        "rxn15121_c0",
    ]
    # Make a copy of the original model to modify
    amac_model_strict_atp = amac_model.copy()
    # Loop through and set all lower bounds to 0
    for rxn_id in atp_consuming_reactions:
        if (rxn_id in amac_model_strict_atp.reactions) & (
            amac_model_strict_atp.metabolites.cpd00002_c0
            in amac_model_strict_atp.reactions.get_by_id(rxn_id).reactants
        ):
            amac_model_strict_atp.reactions.get_by_id(rxn_id).lower_bound = 0
        else:
            print(f"Warning:{rxn_id} not found or ATP is not a reactant.")

    # --- Even stricter: make nucleotide balancing reactions irreversible ---
    # List nucleotide balancing reactions to make irreversible (consumption only)
    nucleotide_balancing_reactions = [
        "rxn00097_c0",
        "rxn00409_c0",
        "rxn00515_c0",
        "rxn00839_c0",
        "rxn01353_c0",
        "rxn01673_c0",
        "rxn01678_c0",
    ]
    # Make a copy of the strict ATP model to modify
    amac_model_strict_nucleotide_balancing = amac_model_strict_atp.copy()
    # Loop through and set all lower bounds to 0
    for rxn_id in nucleotide_balancing_reactions:
        if (rxn_id in amac_model_strict_nucleotide_balancing.reactions) & (
            amac_model_strict_nucleotide_balancing.metabolites.cpd00002_c0
            in amac_model_strict_nucleotide_balancing.reactions.get_by_id(
                rxn_id
            ).reactants
        ):
            amac_model_strict_nucleotide_balancing.reactions.get_by_id(
                rxn_id
            ).lower_bound = 0
        else:
            print(f"Warning:{rxn_id} not found or ATP is not a reactant.")

    # --- The (less) struct ATP production + lumped reactions ---
    amac_model_strict_atp_w_lumped_rxns = amac_model_strict_atp.copy()
    # Remove (not just knockout) the original ICDH and AKGDH reactions
    amac_model_strict_atp_w_lumped_rxns.remove_reactions(
        ["rxn01387_c0", "rxn00199_c0", "rxn00441_c0", "rxn02376_c0", "rxn01872_c0"]
    )
    # Add the lumped reactions
    amac_model_strict_atp_w_lumped_rxns.add_reactions(
        [lumped_icdh_reaction, lumped_akgdh_rxn]
    )
    # --- Make a dictionary of the models to test ---
    amac_models_to_test = {
        "Original": amac_model,
        "With_Lumped_Reactions": amac_w_lumped_rxns,
        "Strict_ATP_Production": amac_model_strict_atp,
        "Strict_ATP_Production_Nucleotide_Balancing": amac_model_strict_nucleotide_balancing,
        "Strict_ATP_Production_With_Lumped_Reactions": amac_model_strict_atp_w_lumped_rxns,
    }

    ####################################################################
    # GROWTH SIMULATIONS (OPTIMIZE FOR BIOMASS)
    ####################################################################
    for model_name, model_obj in amac_models_to_test.items():
        amac_growth_df = run_flux_simulations(
            model=model_obj,
            model_name=model_name,
            basal_media=amac_basal_media,
            c_sources=c_sources,
            n_sources=n_sources,
            id_type="ModelSEED",
            objectives=["bio1_biomass"],
            fluxes_to_record=modelseed_tca_rxns,
        )
        generate_escher_maps(
            df=amac_growth_df,
            model=model_obj,
            map_path=amac_map_path,
            file_prefix=f"amac_{model_name}",
        )

    ecoli_growth_df = run_flux_simulations(
        model=ecoli_model,
        model_name="Ecoli",
        basal_media=ecoli_basal_media,
        c_sources=c_sources,
        n_sources=n_sources,
        id_type="BiGG",
        objectives=["BIOMASS_Ec_iJO1366_core_53p95M"],
        fluxes_to_record=bigg_tca_rxns,
    )
    generate_escher_maps(
        df=ecoli_growth_df,
        model=ecoli_model,
        map_path=ecoli_map_path,
        file_prefix="ecoli",
    )

    ####################################################################
    # BLOCKAGE SIMULATIONS (OPTIMIZE FOR EACH TCA REACTION)
    ####################################################################
    for model_name, model_obj in amac_models_to_test.items():
        amac_blockage_df = run_flux_simulations(
            model=model_obj,
            model_name=model_name,
            basal_media=amac_basal_media,
            c_sources=c_sources,
            n_sources=n_sources,
            id_type="ModelSEED",
            objectives=modelseed_tca_rxns,
            fluxes_to_record=[],
        )
        generate_escher_maps(
            df=amac_blockage_df,
            model=model_obj,
            map_path=amac_map_path,
            file_prefix=f"amac_{model_name}",
            file_suffix="_blocked_reactions",
        )

    ecoli_blockage_df = run_flux_simulations(
        model=ecoli_model,
        model_name="Ecoli",
        basal_media=ecoli_basal_media,
        c_sources=c_sources,
        n_sources=n_sources,
        id_type="BiGG",
        objectives=bigg_tca_rxns,
        fluxes_to_record=[],
    )
    generate_escher_maps(
        df=ecoli_blockage_df,
        model=ecoli_model,
        map_path=ecoli_map_path,
        file_prefix="ecoli",
        file_suffix="_blocked_reactions",
    )


########################################################################
# HELPER FUNCTIONS
########################################################################


def define_media_and_reactions():
    """Returns all media and reaction list definitions."""
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
        "EX_co2_e": 1000,
        "EX_cobalt2_e": 1000,
        "EX_h_e": 1000,
        "EX_h2o_e": 1000,
        "EX_k_e": 1000,
        "EX_cu2_e": 1000,
        "EX_mg2_e": 1000,
        "EX_mn2_e": 1000,
        "EX_mobd_e": 1000,
        "EX_na1_e": 1000,
        "EX_ca2_e": 1000,
        "EX_cbl1_e": 0.01,
        "EX_ni2_e": 1000,
        "EX_o2_e": 1000,
        "EX_cl_e": 1000,
        "EX_pi_e": 1000,
        "EX_zn2_e": 1000,
        "EX_sel_e": 1000,
        "EX_slnt_e": 1000,
        "EX_so4_e": 1000,
        "EX_tungs_e": 1000,
        "EX_fe2_e": 1000,
        "EX_fe3_e": 1000,
    }
    c_sources = {
        "Glucose": {"ModelSEED": "EX_cpd00027_e0", "BiGG": "EX_glc__D_e"},
        "Acetate": {"ModelSEED": "EX_cpd00029_e0", "BiGG": "EX_ac_e"},
    }
    n_sources = {
        "Ammonia": {"ModelSEED": "EX_cpd00013_e0", "BiGG": "EX_nh4_e"},
        "Nitrate": {"ModelSEED": "EX_cpd00209_e0", "BiGG": "EX_no3_e"},
    }
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
    return (
        amac_basal_media,
        ecoli_basal_media,
        c_sources,
        n_sources,
        modelseed_tca_rxns,
        bigg_tca_rxns,
    )


def run_flux_simulations(
    model: cobra.Model,
    model_name: str,
    basal_media: dict,
    c_sources: dict,
    n_sources: dict,
    id_type: str,
    objectives: list,
    fluxes_to_record: list,
) -> pd.DataFrame:
    """Runs FBA for various media conditions and objectives"""
    results_list = []
    c_names = list(c_sources.keys())
    n_names = list(n_sources.keys())

    for c_k in range(1, len(c_names) + 1):
        for c_subset in itertools.combinations(c_names, c_k):
            for n_k in range(1, len(n_names) + 1):
                for n_subset in itertools.combinations(n_names, n_k):
                    c_label = "+".join(c_subset)
                    n_label = "+".join(n_subset)

                    medium = basal_media.copy()
                    for cname in c_subset:
                        medium[c_sources[cname][id_type]] = (
                            10 if cname == "Glucose" else 30
                        )
                    for nname in n_subset:
                        medium[n_sources[nname][id_type]] = 1000

                    m = model.copy()
                    m.medium = medium

                    flux_data = {}
                    for obj_id in objectives:
                        if obj_id not in m.reactions:
                            continue
                        m.objective = obj_id
                        sol = m.optimize()
                        if len(objectives) == 1:
                            flux_data = {
                                rxn: sol.fluxes.get(rxn, 0) for rxn in fluxes_to_record
                            }
                            flux_data[obj_id] = sol.objective_value
                        else:
                            flux_data[obj_id] = sol.objective_value

                    results_list.append(
                        {
                            "Model": model_name,
                            "C_source": c_label,
                            "N_source": n_label,
                            "fluxes": flux_data,
                        }
                    )
    return pd.DataFrame(results_list)


def generate_escher_maps(
    df: pd.DataFrame,
    model: cobra.Model,
    map_path: str,
    file_prefix: str,
    file_suffix: str = "",
):
    """Generates and saves Escher maps from a DataFrame of flux results."""
    for _, row in df.iterrows():
        # Extract the information from the dataframe
        c_label = row["C_source"]
        n_label = row["N_source"]
        flux_data = row["fluxes"]

        # Make a specific file name
        filename = f"{file_prefix}_{c_label}+{n_label}{file_suffix}"

        # Build the Escher Map
        builder = escher.Builder(
            model=model, map_json=map_path, reaction_data=flux_data
        )

        # Save as html
        html_path = os.path.join(HTML_DIR, filename + ".html")
        builder.save_html(html_path)

        # Open in a headless browser, inline computed stules, get SVG
        with sync_playwright() as p:
            browser = p.chromium.launch(headless=True)
            page = browser.new_page(viewport={"width": 1200, "height": 800})
            page.goto(f"file://{html_path.resolve()}", wait_until="domcontentloaded")
            page.wait_for_selector("svg", timeout=10000)

            # Inline computed styles into SVG before serializing
            inline_svg = page.evaluate(
                """() => {
                function copyComputedStyles(svg) {
                    const all = svg.querySelectorAll('*');
                    for (const el of all) {
                        const cs = window.getComputedStyle(el);
                        let s = '';
                        for (let i = 0; i < cs.length; i++) {
                            const key = cs[i];
                            const val = cs.getPropertyValue(key);
                            // Skip some properties if you want smaller output
                            s += key + ':' + val + ';';
                        }
                        el.setAttribute('style', s);
                    }
                }
                const svg = document.querySelector('svg');
                const clone = svg.cloneNode(true);
                copyComputedStyles(clone);
                return new XMLSerializer().serializeToString(clone);
            }"""
            )
            browser.close()

        # Save SVG
        with open(
            os.path.join(SVG_DIR, filename + ".svg"),
            "w",
            encoding="utf-8",
        ) as f:
            f.write(inline_svg)


# Run the main function
if __name__ == "__main__":
    main()
