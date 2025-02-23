import json
import os

import pandas

FILE_DIR = os.path.dirname(os.path.abspath(__file__))

MODELSEED_DIR = "/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase"

# Load the ModelSEED compound database as a dataframe
compounds = pandas.read_csv(
    os.path.join(MODELSEED_DIR, "Biochemistry", "compounds.tsv"),
    sep="\t",
    index_col=0,
)

# Load the ModelSEED Gram Negative biomass template as a DataFrame
biomass = pandas.read_csv(
    os.path.join(MODELSEED_DIR, "Templates", "GramNegative", "BiomassCompounds.tsv"),
    sep="\t",
    index_col=0,
)

# Make the 'id' column the index
biomass.index = biomass["id"]

# Add the compound name, formaula, mass, and charge to the biomass DataFrame
biomass = biomass.join(compounds[["name", "formula", "mass", "charge"]])

# Save the biomass DataFrame as a new TSV file
biomass.to_csv(os.path.join(FILE_DIR, "gram_negative_biomass_components.tsv"), sep="\t")
