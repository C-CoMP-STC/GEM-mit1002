#!/Users/helenscott/opt/miniconda3/envs/escher/bin/python3
"""
This script generates escher maps with fluxes while growing, and while
optimizing for the TCA cycle reactions individually, for the mit1002 model,
the model with some interventions to test, and the E. coli iJO1366 model.
"""
import itertools
import os
import threading
import time
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Optional

import cobra
import pandas as pd
from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
from playwright.sync_api import sync_playwright

import escher

# Define constants and paths at the module level
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(FILE_DIR)
RESULTS_DIR = os.path.join(FILE_DIR, "results")
ESCHER_PLOT_DIR = os.path.join(RESULTS_DIR, "escher_plots")
HTML_DIR = os.path.join(ESCHER_PLOT_DIR, "html")
SVG_DIR = os.path.join(ESCHER_PLOT_DIR, "svg")
PNG_DIR = os.path.join(ESCHER_PLOT_DIR, "png") # Added for PNG output
E_COLI_MODEL_PATH = "/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/ecoli"


def main():
    """
    Main function to run all simulations and generate Escher maps.
    """
    ####################################################################
    # SET UP
    ####################################################################
    # Create directory for escher plots if it doesn't exist
    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(ESCHER_PLOT_DIR, exist_ok=True)
    os.makedirs(HTML_DIR, exist_ok=True)
    os.makedirs(SVG_DIR, exist_ok=True)
    os.makedirs(PNG_DIR, exist_ok=True) # Added for PNG output

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

    # Slightly less strict: Stil allow rxn01517 to be reversible
    # List reactions involving ATP to make irreversible (consumption only)
    atp_consuming_reactions_except_dUMP = [
        "rxn00077_c0",
        "rxn00104_c0",
        "rxn00239_c0",
        "rxn00364_c0",
        "rxn00379_c0",
        "rxn01219_c0",
        "rxn01509_c0",
        "rxn02314_c0",
        "rxn08762_c0",
        "rxn15121_c0",
    ]
    # Make a copy of the original model to modify
    amac_model_strict_atp_except_dUMP = amac_model.copy()
    # Loop through and set all lower bounds to 0
    for rxn_id in atp_consuming_reactions_except_dUMP:
        if (rxn_id in amac_model_strict_atp_except_dUMP.reactions) & (
            amac_model_strict_atp_except_dUMP.metabolites.cpd00002_c0
            in amac_model_strict_atp_except_dUMP.reactions.get_by_id(rxn_id).reactants
        ):
            amac_model_strict_atp_except_dUMP.reactions.get_by_id(rxn_id).lower_bound = 0
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
        "Strict_ATP_Production_Except_dUMP": amac_model_strict_atp_except_dUMP,
        "Strict_ATP_Production_Nucleotide_Balancing": amac_model_strict_nucleotide_balancing,
        "Strict_ATP_Production_With_Lumped_Reactions": amac_model_strict_atp_w_lumped_rxns,
    }

    ####################################################################
    # GROWTH SIMULATIONS (OPTIMIZE FOR BIOMASS)
    ####################################################################
    # For all the different amac models
    for model_name, model_obj in amac_models_to_test.items():
        amac_growth_df = run_flux_simulations(
            model=model_obj,
            model_name=model_name,
            basal_media=amac_basal_media,
            c_sources=c_sources,
            n_sources=n_sources,
            id_type="ModelSEED",
            objectives={"bio1_biomass": "max"},
            fluxes_to_record=[r.id for r in model_obj.reactions],
        )
        # Add a new column with the growth rate
        amac_growth_df["growth_rate"] = amac_growth_df["fluxes"].apply(lambda x: x.get("bio1_biomass", 0))
        # Save the growth simulation results
        amac_growth_df.to_csv(os.path.join(RESULTS_DIR, "amac_" + model_name + "_growth_fluxes.csv"),
                              index=False)
        # Make and save escher maps
        generate_escher_maps(
            df=amac_growth_df,
            model=model_obj,
            map_path=amac_map_path,
            file_prefix=f"amac_{model_name}",
            html_dir=Path(HTML_DIR),
            svg_dir=Path(SVG_DIR),
            png_dir=Path(PNG_DIR),
        )
    # For the one ecoli model
    ecoli_growth_df = run_flux_simulations(
        model=ecoli_model,
        model_name="Ecoli",
        basal_media=ecoli_basal_media,
        c_sources=c_sources,
        n_sources=n_sources,
        id_type="BiGG",
        objectives={"BIOMASS_Ec_iJO1366_core_53p95M": "max"},
        fluxes_to_record=[r.id for r in ecoli_model.reactions],
    )
    # Add a new column with the growth rate
    ecoli_growth_df["growth_rate"] = ecoli_growth_df["fluxes"].apply(lambda x: x.get("BIOMASS_Ec_iJO1366_core_53p95M", 0))
    # Save the growth simulation results
    ecoli_growth_df.to_csv(os.path.join(RESULTS_DIR, "ecoli_growth_fluxes.csv"),
                           index=False)
    # Make and save escher maps
    generate_escher_maps(
        df=ecoli_growth_df,
        model=ecoli_model,
        map_path=ecoli_map_path,
        file_prefix="ecoli",
        html_dir=Path(HTML_DIR),
        svg_dir=Path(SVG_DIR),
        png_dir=Path(PNG_DIR),
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
            html_dir=Path(HTML_DIR),
            svg_dir=Path(SVG_DIR),
            png_dir=Path(PNG_DIR),
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
        html_dir=Path(HTML_DIR),
        svg_dir=Path(SVG_DIR),
        png_dir=Path(PNG_DIR),
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
    # Define the reaction ID for all reactions in the DCA cycle and the
    # direction of flux for the "forward" TCA Cycle (counter clockwise on my
    # map)
    modelseed_tca_rxns = {
        "rxn00256_c0": 'max',  # H2O + Acetyl-CoA + Oxaloacetate --> CoA + H+ + Citrate
        "rxn00974_c0": 'max',  # Citrate <=> H2O + cis-Aconitate
        "rxn01388_c0": 'min',  # Isocitrate <=> H2O + cis-Aconitate
        "rxn01387_c0": 'max',  # NADP + Isocitrate <=> NADPH + H+ + Oxalosuccinate
        "rxn00199_c0": 'max',  # H+ + Oxalosuccinate --> CO2 + 2-Oxoglutarate
        "rxn00198_c0": 'max',  # NADP+ + Isocitrate --> NADPH + H+ + 2-Oxoglutarate
        "rxn00336_c0": 'max',  # Isocitrate --> Succinate + Glyoxalate
        "rxn00330_c0": 'max',  # H2O + Acetyl-CoA + Glyoxalate --> CoA + H+ + L-Malate
        "rxn00441_c0": 'max',  # 2-Oxoglutarate + TPP + H+ --> CO2 + 3-Carboxy-1-hydroxypropyl-TPP
        "rxn02376_c0": 'max',  # 3-Carboxy-1-hydroxypropyl-TPP + Protein N6-(lipoyl)lysine [c0] <=> TPP + S-Succinyldihydrolipoamide
        "rxn01872_c0": 'min',  # Succinyl-CoA + Protein N6-(dihydrolipoyl)lysine <=> CoA + S-Succinyldihydrolipoamide
        "rxn08094_c0": 'max',  # 2-Oxoglutarate + CoA + NAD+ --> Succinyl-CoA + CO2 + NADH
        "rxn00285_c0": 'min',  # ATP + CoA + Succinate <=> ADP + Phosphate + Succinyl-CoA
        "rxn00288_c0": 'max',  # FAD + Succinate + H+ --> Fumarate + FADH2
        "rxn10126_c0": 'max',  # FADH2 + Ubiquinone-8 --> FAD + H+ + Ubiquinol-8
        "rxn00799_c0": 'min',  # L-Malate <=> H2O + Fumarate
        "rxn00248_c0": 'max',  # NAD + L-Malate <=> NADH + Oxaloacetate + H+
    }
    bigg_tca_rxns = {
        "ACONTa": 'max',  # Citrate <=> Cis-Aconitate + H2O H2O
        "ACONTb": 'max',  # Cis-Aconitate + H2O H2O <=> Isocitrate
        "ICDHyr": 'max',  # Isocitrate + NADP <=> 2-Oxoglutarate + CO2 CO2 + NADPH
        "AKGDH": 'max',  # 2-Oxoglutarate + CoA + NAD --> CO2 CO2 + NADH + Succinyl-CoA
        "SUCOAS": 'min',  # 'ATP + CoA + Succinate <=> ADP + Phosphate + Succinyl-CoA'
        "FRD2": 'min',  # Fumarate + Menaquinol 8 --> Menaquinone 8 + Succinate
        "SUCDi": 'max',  # Ubiquinone-8 + Succinate --> Fumarate + Ubiquinol-8
        "FUM": 'max',  # Fumarate + H2O H2O <=> L-Malate
        "MOX": 'max',  # L-Malate + O2 O2 <=> Hydrogen peroxide + Oxaloacetate
        "MDH2": 'max',  # L-Malate + Ubiquinone-8 --> Oxaloacetate + Ubiquinol-8
        "MDH": 'max',  # L-Malate + NAD <=> H+ + NADH + Oxaloacetate
        "MDH3": 'max',  # L-Malate + Menaquinone 8 --> Menaquinol 8 + Oxaloacetate
        "CS": 'max',  # Acetyl-CoA + H2O H2O + Oxaloacetate --> Citrate + Coenzyme A + H+
        "CITL": 'min',  # Citrate --> Acetate + Oxaloacetate
        "MALS": 'max',  # Acetyl-CoA + Glyoxylate + H2O H2O --> Coenzyme A + H+ + L-Malate
        "ICL": 'max',  # Isocitrate --> Glyoxylate + Succinate'
    }
    return (
        amac_basal_media,
        ecoli_basal_media,
        c_sources,
        n_sources,
        modelseed_tca_rxns,
        bigg_tca_rxns,
    )


# Quiet HTTP handler so the server doesn't spam stdout
class _QuietHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        # Use the 'directory' keyword argument for Python 3.7+
        super().__init__(*args, directory=kwargs.pop("directory"), **kwargs)

    def log_message(self, format, *args):
        pass


def _start_http_server(directory: Path, port: int = 0):
    """Starts an HTTP server in a thread, serving the specified directory."""
    # This handler is created with the directory to serve, avoiding os.chdir
    handler = lambda *args, **kwargs: _QuietHandler(*args, directory=str(directory), **kwargs)
    server = ThreadingHTTPServer(("127.0.0.1", port), handler)
    actual_port = server.server_address[1]
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    # Return only the three values expected by the calling function
    return server, thread, actual_port


def run_flux_simulations(
    model: cobra.Model,
    model_name: str,
    basal_media: dict,
    c_sources: dict,
    n_sources: dict,
    id_type: str,
    objectives: dict,
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
                    for obj_id, obj_dir in objectives.items():
                        if obj_id not in m.reactions:
                            continue
                        m.objective = obj_id
                        m.objective_direction = obj_dir
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
    html_dir: Optional[Path] = None,
    svg_dir: Optional[Path] = None,
    png_dir: Optional[Path] = None,
    min_svg_inner_length: int = 700,
    page_viewport: dict = {"width": 1200, "height": 800},
) -> None:
    """
    Generates and saves Escher maps from a DataFrame of flux results.
    Saves standalone SVGs and PNGs to disk.
    """
    # Set default directories if not provided
    html_dir = Path(html_dir) if html_dir else Path("html")
    svg_dir = Path(svg_dir) if svg_dir else Path("svg")
    png_dir = Path(png_dir) if png_dir else Path("png")

    html_dir.mkdir(parents=True, exist_ok=True)
    svg_dir.mkdir(parents=True, exist_ok=True)
    png_dir.mkdir(parents=True, exist_ok=True)

    # Start HTTP server serving the html_dir
    server, thread, port = _start_http_server(html_dir, port=0)
    base_url = f"http://127.0.0.1:{port}"
    print(f"Serving {html_dir.resolve()} at {base_url}/")

    try:
        with sync_playwright() as p:
            browser = p.chromium.launch(headless=True)
            page = browser.new_page(viewport=page_viewport)
            page.on("console", lambda msg: print(f"BROWSER: {msg.text}"))

            for _, row in df.iterrows():
                c_label = row["C_source"]
                n_label = row["N_source"]
                flux_data = row["fluxes"]

                filename = f"{file_prefix}_{c_label}+{n_label}{file_suffix}"
                html_filename = filename + ".html"
                svg_path = svg_dir / (filename + ".svg")
                png_path = png_dir / (filename + ".png")

                builder = escher.Builder(
                    model=model, map_json=map_path, reaction_data=flux_data
                )
                # Save using the full path now that we are not changing directory
                builder.save_html(str(html_dir / html_filename))

                url = f"{base_url}/{html_filename}"
                print("Loading", url)

                try:
                    page.goto(url, wait_until="domcontentloaded", timeout=15000)
                    page.wait_for_selector("svg.escher-svg", timeout=15000)
                    page.wait_for_function(
                        f"() => {{ const s = document.querySelector('svg.escher-svg'); return s && s.innerHTML.length > {min_svg_inner_length}; }}",
                        timeout=20000,
                    )
                except PlaywrightTimeoutError as e:
                    print(f"Warning: Timeout loading page or SVG for {filename}. Error: {e}. Skipping image generation for this file.")
                    continue

                # --- Generate SVG ---
                inline_svg = page.evaluate(
                    """() => {
                    function inlineStyles(svg) {
                        for (const el of svg.querySelectorAll('*')) {
                            const cs = window.getComputedStyle(el);
                            if (cs.length > 0) el.setAttribute('style', cs.cssText);
                        }
                    }
                    const svg = document.querySelector('svg.escher-svg');
                    if (!svg) return '';
                    const clone = svg.cloneNode(true);
                    inlineStyles(clone);
                    if (!clone.hasAttribute('xmlns')) clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
                    return new XMLSerializer().serializeToString(clone);
                }"""
                )
                if inline_svg:
                    svg_path.write_text(inline_svg, encoding="utf-8")
                    print("Wrote", svg_path)
                else:
                    print(f"Warning: extracted SVG for {filename} was empty.")

                # --- Generate PNG ---
                try:
                    el = page.query_selector('svg.escher-svg')
                    if el:
                        el.screenshot(path=str(png_path))
                        print("Wrote", png_path)
                except Exception as e:
                    print(f"PNG screenshot failed for {filename}: {e}")

            browser.close()
    finally:
        server.shutdown()
        print("HTTP server stopped.")


# Run the main function
if __name__ == "__main__":
    main()
