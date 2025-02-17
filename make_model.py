#!~/opt/miniconda3/envs/modelseedpy-dev/bin/python
"""
Script to build a ModelSEED model from a genome, add biomass and gap‐filled reactions,
and save intermediate SBML files.
"""

import ast
import json
import os
import pickle

import cobra
import pandas as pd
from modelseedpy import FBAHelper, KBaseMediaPkg, MSBuilder, MSGenome, RastClient
from modelseedpy.core.msbuilder import build_biomass, core_atp
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.helpers import get_template

# Define directories
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
TESTFILE_DIR = os.path.join(FILE_DIR, "test", "test_files")


# =============================================================================
#   Load Genome and Annotations
# =============================================================================
# Load the genome (from Michelle's anvi'o gene calls)
genome = MSGenome.from_fasta(
    "genome/Michelle's 4106 gene calls/MIT1002_anvio_prot_seqs.fa", split=" "
)


def parse_gff(gff_path):
    """
    Parses a GFF file and returns a dictionary mapping gene IDs to annotations.

    Args:
        gff_path (str): Path to the GFF file

    Returns:
        dict: Mapping of gene ID to its annotation dictionary.
    """
    gene_annotations = {}
    with open(gff_path, "r") as gff_file:
        for line in gff_file:
            if line.startswith("#"):
                continue  # Skip comments
            columns = line.strip().split("\t")
            if len(columns) < 9:
                continue  # Skip malformed lines

            feature_type = columns[2]
            attributes = columns[8]
            gene_id = None
            for attr in attributes.split(";"):
                if attr.startswith("ID=") or attr.startswith("locus_tag="):
                    gene_id = attr.split("=")[1]
                    break
            if gene_id:
                gene_annotations[gene_id] = {
                    "feature_type": feature_type,
                    "info": attributes,
                }
    return gene_annotations


# Load GFF annotations and attach them to genome features
gff_annotations = parse_gff(
    "genome/Michelle's 4106 gene calls/2738541267_genecalls.gff"
)
for feature in genome.features:
    if feature.id in gff_annotations:
        feature.annotation = gff_annotations[feature.id]

# Inspect an example feature
print(genome.features[0].id, genome.features[0].annotation)

# Annotate the genome using RAST
rast = RastClient()
rast.annotate_genome(genome)


# =============================================================================
#   Build the Model and Add Biomass Reaction
# =============================================================================
# Use the gram-negative template (see FIXME below if you later want to classify the genome)
template = get_template("template_gram_neg")
cobra_template = MSTemplateBuilder.from_dict(template, None).build()

# Build the initial model from the genome and template
model_builder = MSBuilder(genome, cobra_template)
base_model = model_builder.build(
    "model",
    "0",
    allow_all_non_grp_reactions=True,
    annotate_with_rast=False,
)

# Add a biomass reaction (using biomass ID "bio2")
# FIXME: Check that this biomass formulation is correct for your application.
base_model.add_reactions(
    [build_biomass("bio2", base_model, cobra_template, core_atp, "0")]
)

# Save the base model
cobra.io.write_sbml_model(base_model, "modelseedpy_model_01.xml")


# =============================================================================
#   Add Michelle's Reactions from the ModelSEED Database
# =============================================================================
# Load ModelSEED reaction and compound databases
rxn_db = json.load(
    open(
        "/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase/Biochemistry/reactions.json"
    )
)
met_db = json.load(
    open(
        "/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase/Biochemistry/compounds.json"
    )
)

# Get the list of reaction IDs from the template (remove the trailing compartment tag)
template_rxn_ids = [r["id"][:-2] for r in template["reactions"]]

# Subset the ModelSEED reaction DB to non-obsolete reactions present in the template
template_rxn_db = {
    rxn["id"]: rxn
    for rxn in rxn_db
    if not rxn["is_obsolete"] and rxn["id"] in template_rxn_ids
}

# Load Michelle's reactions (with ModelSEED IDs) from CSV
michelle_rxns = pd.read_csv(
    os.path.join("Pangenome from Michelle", "database_w_MNX_SEED.csv"), header=0
)

# Filter for inferred reactions with a valid ModelSEED ID and convert the IDs from string to list
rxns_to_add = michelle_rxns[
    (michelle_rxns["Inferred Presence"] == 1)
    & (michelle_rxns["ModelSEED ID"].notnull())
].copy()
rxns_to_add["ModelSEED ID"] = rxns_to_add["ModelSEED ID"].apply(ast.literal_eval)

# Combine the lists of reaction IDs and get unique values, then subset to those in the template DB
lists_of_rxn_ids = rxns_to_add["ModelSEED ID"].tolist()
rxn_ids = list(
    {
        x
        for item in lists_of_rxn_ids
        for x in (item if isinstance(item, list) else [item])
    }
)
rxn_ids = [rxn_id for rxn_id in rxn_ids if rxn_id in template_rxn_db]


def create_cobra_reaction(model, modelseed_db, rxn_id):
    """
    Create a COBRApy Reaction object from a ModelSEED database entry.

    Parameters:
        model (cobra.Model): Model to which the reaction and new metabolites will be added.
        modelseed_db (dict): Dictionary of ModelSEED reaction entries.
        rxn_id (str): Reaction ID to add.

    Returns:
        cobra.Reaction: The created reaction.
    """
    rxn = modelseed_db[rxn_id]
    reaction_id = rxn["id"] + "_c0"  # Assume cytosolic compartment
    reaction = cobra.Reaction(reaction_id, name=rxn["name"])

    # Set reaction bounds based on reversibility
    if rxn["reversibility"] == ">":
        reaction.lower_bound = 0
        reaction.upper_bound = 1000
    elif rxn["reversibility"] == "<":
        reaction.lower_bound = -1000
        reaction.upper_bound = 0
    else:  # '=' or any unexpected value
        reaction.lower_bound = -1000
        reaction.upper_bound = 1000

    # Add metabolites to the reaction
    for met_str in rxn["stoichiometry"].split(";"):
        parts = met_str.split(":")
        if len(parts) < 3:
            continue  # Skip malformed entries

        coeff = float(parts[0])
        met_id_base = parts[1]
        comp_code = parts[2]

        # Map compartment code to tag
        if comp_code == "0":
            compartment = "c0"
        elif comp_code == "2":
            compartment = "e0"
        else:
            compartment = comp_code

        met_id = f"{met_id_base}_{compartment}"
        met_name = parts[4].strip('"') if len(parts) >= 5 else met_id

        # Add metabolite if not present
        if met_id not in model.metabolites:
            met_obj = cobra.Metabolite(
                id=met_id, name=met_name, compartment=compartment
            )
            model.add_metabolites([met_obj])
        else:
            met_obj = model.metabolites.get_by_id(met_id)

        reaction.add_metabolites({met_obj: coeff})

    model.add_reactions([reaction])
    return reaction


# Add each selected reaction to the base model
for rxn_id in rxn_ids:
    create_cobra_reaction(base_model, template_rxn_db, rxn_id)

# Save the updated model
cobra.io.write_sbml_model(base_model, "modelseedpy_model_02.xml")


# =============================================================================
#   Gapfill and Annotate Biomass Components
# =============================================================================
# Load media definitions
with open(os.path.join(TESTFILE_DIR, "media", "media_definitions.pkl"), "rb") as f:
    media_definitions = pickle.load(f)

# Add glucose to the mbm media
media_definitions["mbm_media"]["EX_cpd00027_e0"] = 10


def gapfill_and_annotate_biomass_components(model, template, media, biomass_rxn_id):
    """
    For each biomass component in the biomass reaction, run gap filling by adding a temporary sink
    reaction for that component. The final model is the union of all gap-filled reactions, and each
    reaction is annotated with a dictionary mapping biomass component IDs to a binary flag
    (1 = added by gap filling, 0 = not added).

    Parameters:
        model (cobra.Model): The original COBRApy model.
        template: Template (e.g., from KBase) for gap filling.
        media: Media object for gap filling.
        biomass_rxn_id (str): ID of the biomass reaction.

    Returns:
        cobra.Model: The final gap-filled model.
    """
    # Start with a copy of the original model and add sink reactions for all metabolites.
    base_model = model.copy()
    for metabolite in model.metabolites:
        base_model.add_boundary(metabolite, type="sink", lb=0)

    try:
        biomass_rxn = base_model.reactions.get_by_id(biomass_rxn_id)
    except KeyError:
        raise ValueError(f"Biomass reaction {biomass_rxn_id} not found in the model.")

    biomass_components = list(biomass_rxn.metabolites.keys())
    gapfill_results = {}

    # Create a union model to collect all gap-filled reactions.
    final_model = base_model.copy()

    # Initialize each reaction’s gap-fill annotation.
    for rxn in final_model.reactions:
        rxn.annotation.setdefault("gapfill_results", {})
        for met in biomass_components:
            rxn.annotation["gapfill_results"].setdefault(met.id, 0)

    # Process each biomass component individually.
    for met in biomass_components:
        temp_model = final_model.copy()
        sink_id = f"SK_{met.id}"
        temp_model.objective = sink_id

        # Run gap filling (MSBuilder.gapfill_model returns a gap-filled model).
        gapfilled_model = MSBuilder.gapfill_model(temp_model, sink_id, template, media)

        # Determine added reactions by comparing reaction IDs.
        temp_rxn_ids = {r.id for r in temp_model.reactions}
        added_rxn_ids = [
            r.id for r in gapfilled_model.reactions if r.id not in temp_rxn_ids
        ]
        gapfill_results[met.id] = added_rxn_ids

        # Add new gap-filled reactions to the final model and annotate them.
        for rxn_id in added_rxn_ids:
            if rxn_id not in final_model.reactions:
                gap_rxn = gapfilled_model.reactions.get_by_id(rxn_id)
                gap_rxn.annotation.setdefault("gapfill_results", {})
                for bm in biomass_components:
                    gap_rxn.annotation["gapfill_results"].setdefault(bm.id, 0)
                final_model.add_reactions([gap_rxn.copy()])

        # Mark reactions added in this gap-fill run.
        for rxn in final_model.reactions:
            if rxn.id in added_rxn_ids:
                rxn.annotation["gapfill_results"][met.id] = 1

    # Optionally remove temporary sink reactions.
    sk_rxn_ids = [rxn.id for rxn in final_model.reactions if rxn.id.startswith("SK_")]
    if sk_rxn_ids:
        final_model.remove_reactions(sk_rxn_ids, remove_orphans=True)

    return final_model


# # Run gap filling for each biomass component in the biomass reaction "bio1"
# final_model = gapfill_and_annotate_biomass_components(
#     base_model, template, media_definitions["mbm_media"], "bio1"
# )

# for rxn in final_model.reactions:
#     if "gapfill_results" in rxn.annotation:
#         # Convert the gapfill_results dictionary into a JSON string
#         rxn.annotation["gapfill_results"] = json.dumps(
#             rxn.annotation["gapfill_results"]
#         )

# # Save the final gap-filled model
# cobra.io.write_sbml_model(final_model, "modelseedpy_model_03.xml")

# =============================================================================
#   Gapfill for growth (really just testing the gapfilling)
# =============================================================================
gapfilled_model = MSBuilder.gapfill_model(
    base_model, "bio1", template, media_definitions["mbm_media"]
)

# Save the final gap-filled model
cobra.io.write_sbml_model(gapfilled_model, "modelseedpy_model_04.xml")
