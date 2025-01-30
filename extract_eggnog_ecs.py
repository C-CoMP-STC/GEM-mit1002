import os

import pandas as pd

# Load the eggnong annotations output tsv as a dataframe
eggnog_annotations = pd.read_csv(
    os.path.join(os.path.dirname(__file__), "genome", "eggnog_output.emapper.annotations"), sep="\t", comment="#"
)

eggnog_annotations