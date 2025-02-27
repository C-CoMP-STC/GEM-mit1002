import os

import pandas

FILE_DIR = os.path.dirname(os.path.abspath(__file__))

# Load the eggnog TSV file
eggnog_path = os.path.join(FILE_DIR, 'eggnog_output.emapper.annotations')
eggnog_df = pandas.read_csv(eggnog_path, sep='\t', comment='#', header=None)
eggnog_df.columns = [
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

# Add "anvio_prot_seqs" before every entry in the first "query" column
eggnog_df["query"] = eggnog_df["query"].apply(lambda x: "anvio_prot_seqs_" + str(x))

# Save to a new file
output_path = os.path.join(FILE_DIR, 'clean_eggnog_output.emapper.annotations')
eggnog_df.to_csv(output_path, sep='\t', index=False)
