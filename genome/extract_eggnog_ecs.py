import os

import pandas as pd

# Load the eggnong annotations output tsv as a dataframe
eggnog_annotations = pd.read_csv(
    os.path.join(
        os.path.dirname(__file__), "genome", "eggnog_output.emapper.annotations"
    ),
    sep="\t",
    comment="#",
)
eggnog_annotations.columns = [
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "max_annot_lvl",
    "COG_category",
    "Description",
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
]

# Get the EC column
eggnog_annotations["EC"]

# Print the EC numbers (to share with Daniel)
print(list(eggnog_annotations["EC"].unique()))

# Check the number of bigg reactions in the dataset
eggnog_annotations["BiGG_Reaction"].nunique()
# 421 unique BiGG reactions
# That seems like it should be enough...
