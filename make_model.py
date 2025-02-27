#!/projectnb/cometsfba/hscott/GEM-repos/GEM-mit1002/env/bin/python
"""
Script to build a draft ModelSEED model from a genome annotated with RAST
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

def main():

    # =============================================================================
    #   Load Genome and Annotations
    # =============================================================================
    # Load the genome (from Michelle's anvi'o gene calls)
    genome = MSGenome.from_fasta(
        os.path.join(FILE_DIR, "genome", "Michelle's 4106 gene calls", "MIT1002_anvio_prot_seqs.fa"),
        split=" "
    )

    # Annotate the genome using RAST
    rast = RastClient()
    rast.annotate_genome(genome)

    # =============================================================================
    #   Build the Model and Add Biomass Reaction
    # =============================================================================
    # Use the gram-negative template
    template = get_template("template_gram_neg")
    cobra_template = MSTemplateBuilder.from_dict(template, None).build()

    # Build the initial model from the genome and template
    model_builder = MSBuilder(genome, cobra_template)
    base_model = model_builder.build(
        "iHS4106",  # Model ID
        allow_all_non_grp_reactions=True,
        annotate_with_rast=False,
    )

    # Save the base model
    cobra.io.write_sbml_model(base_model, os.path.join(FILE_DIR, "model.xml"))


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


# Define my own classed to match the KBaseMediaPkg class
class MediaCompound:
    def __init__(self, compound_id, maxFlux, minFlux, concentration=0.001):
        self.id = compound_id
        self.maxFlux = maxFlux
        self.minFlux = minFlux
        self.concentration = concentration


class CustomMedia:
    def __init__(self, name, compounds):
        self.name = name
        self.mediacompounds = compounds


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

    # Process each biomass component individually.
    for met in biomass_components:
        temp_model = base_model.copy()
        sink_id = f"SK_{met.id}"
        temp_model.objective = sink_id

        # Run gap filling (MSBuilder.gapfill_model returns a gap-filled model).
        gapfilled_model = MSBuilder.gapfill_model(temp_model, sink_id, template, media)

        # TODO: Should I remove the sink reactions before saving?

        # Save the gapfilled model
        cobra.io.write_sbml_model(
            gapfilled_model,
            f"modelseedpy_gapfill_per_biomass_cmpt/gapfilled_{met.id}.xml",
        )


if __name__ == "__main__":
    main()