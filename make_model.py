import ast
import json
import os
import pickle
import sys

import cobra
import pandas as pd
from modelseedpy import FBAHelper, KBaseMediaPkg, MSBuilder, MSGenome, RastClient
from modelseedpy.core.msbuilder import build_biomass, core_atp
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.helpers import get_classifier, get_template

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
TESTFILE_DIR = os.path.join(FILE_DIR, "test", "test_files")

# Load the genome (From Michelle's anvi'o gene calls)
genome = MSGenome.from_fasta(
    "genome/Michelle's 4106 gene calls/MIT1002_anvio_prot_seqs.fa", split=" "
)

# Add Michelle's anotations
# Right now, this can't be used in the model building process, but I still want to load them
########################################################################


def parse_gff(gff_path):
    """Parses a GFF file and returns a dictionary mapping gene IDs to annotations.
    Written by ChatGPT

    Args:
        gff_path (str): Path to the GFF file

    Returns:
        dict: A dictionary mapping gene IDs to annotations
    """
    gene_annotations = {}

    with open(gff_path, "r") as gff_file:
        for line in gff_file:
            if line.startswith("#"):
                continue  # Skip comments

            columns = line.strip().split("\t")
            if len(columns) < 9:
                continue  # Skip malformed lines

            # FIXME: Remove magic numbers
            feature_type = columns[2]  # e.g., "CDS", "gene", "mRNA"
            attributes = columns[8]  # This contains gene IDs and other info

            # Extract gene ID (format depends on annotation system)
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


gff_annotations = parse_gff(
    "genome/Michelle's 4106 gene calls/2738541267_genecalls.gff"
)

for feature in genome.features:
    if feature.id in gff_annotations:
        # FIXME: Should I be putting this in ontology_term instead
        feature.annotation = gff_annotations[feature.id]  # Add GFF metadata

# insepct an example feature
print(genome.features[0].id, genome.features[0].annotation)

########################################################################

# Annotate the genome
# FIXME: Do I need this? Can I use the exact same annotations as Michelle?
rast = RastClient()
rast.annotate_genome(genome)

# Export the annotated genome to be able to compare it to Michelle's?

# # Classify the genome
# # Assuming it is classifying the genome as gram positive or negative
# # N would mean gram negative
# genome_classifier = get_classifier("knn_filter")
# genome_class = genome_classifier.classify(genome)

# # Get the template model
# if genome_class == "G":
#     # Not actually sure if G is gram positive, but that is what copilot says
#     template = get_template("template_gram_pos")
# elif genome_class == "N":
#     template = get_template("template_gram_neg")
# else:
#     # Throw an error if the classifier returns something unexpected
#     raise ValueError(f"Unexpected genome class: {genome_class}")

# Getting a weird error within "get_classifier" that I don't understand
# Traceback (most recent call last):
#   File "/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/mit1002-model/make_model.py", line 26, in <module>
#     genome_classifier = get_classifier("knn_filter")
#   File "/Users/helenscott/opt/miniconda3/envs/gem-reconstruction/lib/python3.8/site-packages/modelseedpy/helpers.py", line 43, in get_classifier
#     model_filter = pickle.load(fh)
# ModuleNotFoundError: No module named 'sklearn.neighbors._dist_metrics'
# So I'm just going to define that the template is the gram negative one, since I know that
template = get_template("template_gram_neg")

# Convert to Cobra?
# Not sure what this is doing
cobra_template = MSTemplateBuilder.from_dict(template, None).build()

# Make the model
model_builder = MSBuilder(genome, cobra_template)
base_model = model_builder.build(
    # TODO: Look up what these options are
    "model",
    "0",
    allow_all_non_grp_reactions=True,
    annotate_with_rast=False,
)

# Add the biomass reaction?
# Not sure what this is doing
# FIXME: Check the template is right (in the demo was using core template)
# TODO: Can I make my own biomass reaction? I want to add PHB to it
# TODO: Add another core/minimal biomass reaction
base_model.add_reactions(
    [build_biomass("bio2", base_model, cobra_template, core_atp, "0")]
)

# Save the model
cobra.io.write_sbml_model(base_model, "modelseedpy_model_01.xml")

# Add all of Michelle's reactions
# Need a helper function to buid a COBRA reaction from given a ModelSEED reaction ID and compartment (using the ModelSEED database)
# Need to load the ModelSEED database first
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
# Get all of the reaction IDs in the template
# -2 is a magic number to get rid of the ""_c" at the end of the reaction ID
template_rxn_ids = [r["id"][:-2] for r in template["reactions"]]
# Subset the rxn_db to only include non-obsolete reactions, and reactions in template
# And make it into a dictionary for easier searching
template_rxn_db = {
    rxn["id"]: rxn
    for rxn in rxn_db
    if not rxn["is_obsolete"] and rxn["id"] in template_rxn_ids
}
# Load the database of Michelle's reactions with the ModelSEED IDs
michelle_rxns = pd.read_csv(
    os.path.join("Pangenome from Michelle", "database_w_MNX_SEED.csv"),
    header=0,
)

# Filter the list to only include ones with a "1" for "Inferred Presence" and with something in "ModelSEED ID"
# Make an explicit copy after filtering to avoid SettingWithCopyWarning
rxns_to_add = michelle_rxns[
    (michelle_rxns["Inferred Presence"] == 1)
    & (michelle_rxns["ModelSEED ID"].notnull())
].copy()

# Check if any of the "ModelSEED ID" values are empty
# rxns_to_add["ModelSEED ID"].isnull().sum()

# Convert the "ModelSEED ID" column from a string to a list
rxns_to_add["ModelSEED ID"] = rxns_to_add["ModelSEED ID"].apply(ast.literal_eval)

# Combine all of the values in the "ModelSEED ID" column into a list
lists_of_rxn_ids = rxns_to_add["ModelSEED ID"].tolist()
rxn_ids = [
    x for item in lists_of_rxn_ids for x in (item if isinstance(item, list) else [item])
]
# Get just the unique values
rxn_ids = list(set(rxn_ids))

# Subset the rxn_ids list to only include those in template_rxn_db
rxn_ids = [rxn_id for rxn_id in rxn_ids if rxn_id in template_rxn_db]


def create_cobra_reaction(model, modelseed_db, rxn_id):
    """
    Create a COBRApy Reaction object from a ModelSEED database entry.

    Parameters:
        model (cobra.Model): A COBRApy model to which the reaction and any new metabolites will be added.
        modelseed_db (dict): A dictionary with ModelSEED reaction entries.
        rxn_id (str): The ID of the reaction in the ModelSEED database.

    Returns:
        cobra.Reaction: The constructed reaction.
    """
    rxn = modelseed_db[rxn_id]

    # Create a reaction ID by appending a default compartment tag.
    # (Here we assume the reaction is localized to the cytosol.)
    reaction_id = rxn["id"] + "_c0"
    reaction = cobra.Reaction(reaction_id, name=rxn["name"])

    # --- Handle reversibility ---
    # Mapping the ModelSEED "reversibility" field to bounds:
    #   '>' means irreversible forward: [0, 1000]
    #   '<' means irreversible reverse: [-1000, 0]
    #   '=' means reversible: [-1000, 1000]
    if rxn["reversibility"] == ">":
        reaction.lower_bound = 0
        reaction.upper_bound = 1000
    elif rxn["reversibility"] == "<":
        reaction.lower_bound = -1000
        reaction.upper_bound = 0
    elif rxn["reversibility"] == "=":
        reaction.lower_bound = -1000
        reaction.upper_bound = 1000
    else:
        # Default: reversible reaction
        reaction.lower_bound = -1000
        reaction.upper_bound = 1000

    # --- Add metabolites to the reaction ---
    # The stoichiometry field is a semicolon-delimited string. Each term has the format:
    # "coefficient:met_id:comp_code:...:met_name"
    # Example: '-1:cpd00001:0:0:"H2O"'
    for met_str in rxn["stoichiometry"].split(";"):
        parts = met_str.split(":")
        if len(parts) < 3:
            # Skip malformed entries
            continue

        coeff = float(parts[0])
        met_id_base = parts[1]
        comp_code = parts[2]

        # Determine compartment tag.
        # Here we use '0' to represent the cytosol (c0) and '2' for the extracellular space (e0).
        # You can extend this mapping if necessary.
        if comp_code == "0":
            compartment = "c0"
        elif comp_code == "2":
            compartment = "e0"
        else:
            compartment = comp_code  # Fallback: use the code as given

        met_id = met_id_base + "_" + compartment

        # Use the metabolite name if provided (usually in parts[4]); otherwise, default to met_id.
        if len(parts) >= 5:
            met_name = parts[4].strip('"')
        else:
            met_name = met_id

        # If the metabolite is not already in the model, add it.
        if met_id not in model.metabolites:
            met_obj = cobra.Metabolite(
                id=met_id, name=met_name, compartment=compartment
            )
            model.add_metabolites([met_obj])
        else:
            met_obj = model.metabolites.get_by_id(met_id)

        # Add the metabolite with its stoichiometric coefficient to the reaction.
        reaction.add_metabolites({met_obj: coeff})

    # Add the reaction to the model.
    model.add_reactions([reaction])

    return reaction


# Add all of the reactions to the model
for rxn_id in rxn_ids:
    if rxn_id in template_rxn_db:
        create_cobra_reaction(base_model, template_rxn_db, rxn_id)

# Save the model
cobra.io.write_sbml_model(base_model, "modelseedpy_model_02.xml")

# Load my media
# Load the media definitions
with open(os.path.join(TESTFILE_DIR, "media", "media_definitions.pkl"), "rb") as f:
    media_definitions = pickle.load(f)


def gapfill_and_annotate_biomass_components(model, template, media, biomass_rxn_id):
    """
    For each biomass component in the biomass reaction, run gap filling by
    adding a temporary demand reaction for that component. Then, union all gap filled
    reactions into the final model and update each reaction's annotation with a
    dictionary mapping biomass component IDs to a binary flag (1 = gap filled
    for that component, 0 = not gap filled).

    Parameters:
      model (cobra.Model): The original COBRApy model.
      template: A template (e.g., from KBase) needed by MSBuilder.gapfill_model.
      media: A media object (e.g., from KBaseMedia).
      biomass_rxn_id (str): The ID of the biomass reaction in the model.

    Returns:
      cobra.Model: The final model containing the union of gap filled reactions,
                   with each reaction annotated with its gap fill results.
    """
    # Make a copy of the original model to serve as the base.
    base_model = model.copy()

    # Add sink reactions for all metabolites, but set the lower bound to 0
    # because by default the sink reactions are reversible, and so can be
    # used to import metabolites that are not in the media
    for metabolite in model.metabolites:
        base_model.add_boundary(metabolite, type="sink", lb=0)

    # Get the biomass reaction and the list of biomass components (metabolites).
    try:
        biomass_rxn = base_model.reactions.get_by_id(biomass_rxn_id)
    except KeyError:
        raise ValueError(f"Biomass reaction {biomass_rxn_id} not found in the model.")
    biomass_components = list(biomass_rxn.metabolites.keys())

    # Prepare a dictionary to store gap filling results for each biomass component.
    # Each key will be a biomass metabolite ID and the value will be a list of reaction IDs
    # that were added by gap filling.
    gapfill_results = {}

    # For unioning gap filled reactions we start with the base model.
    final_model = base_model.copy()

    # Initialize each reaction's annotation field for gap fill results.
    # We will add a dictionary under reaction.annotation["gapfill_results"] where keys
    # are biomass component metabolite IDs and values are 0 (by default).
    for rxn in final_model.reactions:
        # If the reaction already has an annotation dict, leave it;
        # otherwise, initialize an empty dict.
        rxn.annotation.setdefault("gapfill_results", {})
        for met in biomass_components:
            # Use the metabolite's ID as key.
            rxn.annotation["gapfill_results"].setdefault(met.id, 0)

    # Process each biomass component individually.
    for met in biomass_components:
        # Create a temporary model that starts from the current union (final_model).
        temp_model = final_model.copy()

        # Use the sink reaction for the current biomass component as the objective
        sink_id = f"SK_{met.id}"
        temp_model.objective = sink_id

        # Run gap filling on the temporary model for this demand reaction.
        # (MSBuilder.gapfill_model returns a gap filled model.)
        gapfilled_model = MSBuilder.gapfill_model(temp_model, sink_id, template, media)

        # Identify which reactions were added by comparing gapfilled_model to temp_model.
        temp_rxn_ids = {r.id for r in temp_model.reactions}
        added_rxn_ids = [
            r.id for r in gapfilled_model.reactions if r.id not in temp_rxn_ids
        ]
        gapfill_results[met.id] = added_rxn_ids

        # Add any gap filled reactions that are new to the union (final_model).
        for rxn_id in added_rxn_ids:
            if rxn_id not in final_model.reactions:
                # Copy the reaction from the gap filled model.
                gap_rxn = gapfilled_model.reactions.get_by_id(rxn_id)
                # Ensure that the gapfill_results annotation exists.
                gap_rxn.annotation.setdefault("gapfill_results", {})
                # Initialize the binary flags for all biomass components.
                for bm in biomass_components:
                    gap_rxn.annotation["gapfill_results"].setdefault(bm.id, 0)
                # Add the reaction to the final model.
                final_model.add_reactions([gap_rxn.copy()])

        # For all reactions in the final model, mark those that were added in this gap fill
        # run (i.e. their IDs appear in added_rxn_ids) with a 1 for this biomass component.
        for rxn in final_model.reactions:
            # (We already initialized each reaction's gapfill_results for this met to 0.)
            if rxn.id in added_rxn_ids:
                rxn.annotation["gapfill_results"][met.id] = 1

    # (Optionally, you could remove the temporary demand reactions from final_model.)
    # For example, to remove all reactions with IDs starting with "DM_":
    sk_rxn_ids = [rxn.id for rxn in final_model.reactions if rxn.id.startswith("SK_")]
    if sk_rxn_ids:
        final_model.remove_reactions(sk_rxn_ids, remove_orphans=True)

    # Now, each reaction in final_model.annotation["gapfill_results"] is a dictionary, for example:
    #   {"cpd00001_c0": 1, "cpd00012_c0": 0, ...}
    # which tells you whether that reaction was added during gap filling for each biomass component.
    return final_model


# =============================================================================
# === Example usage ==========================================================
# =============================================================================
# Assume you already have a COBRApy model, a template (e.g., from KBase),
# and a media object (e.g., from KBaseMedia):
#
#     from modelseedpy import MSBuilder
#     model = ...      # your existing model
#     template = ...   # e.g., kbase_api.get_from_ws("GramNegModelTemplateV3", "NewKBaseModelTemplates")
#     media = ...      # e.g., kbase_api.get_from_ws("Carbon-D-Glucose", "KBaseMedia")
#
# And your biomass reaction ID is "bio1". Then, you can do:
#
final_model = gapfill_and_annotate_biomass_components(
    base_model, template, media_definitions["mbm_media"], "bio1"
)
#
# Now, for any reaction in final_model, you can inspect:
#    print(final_model.reactions.get_by_id("R_12DGR120tipp").annotation["gapfill_results"])
#
# which might print something like:
#    {'cpd00001_c0': 1, 'cpd00012_c0': 0, ...}
#
# This final model (which is the union of all gap filled reactions) can be saved to file,
# and each reaction’s annotation carries the binary ensemble state information.

# Save the final model
cobra.io.write_sbml_model(final_model, "modelseedpy_model_02.xml")
