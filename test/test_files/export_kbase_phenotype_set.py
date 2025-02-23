import os

import pandas as pd

# TODO: Set a variable for the one kbase workspace where all my media are stored
media_workspace = "hgsco:narrative_1726668344695"

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
media = {
    "mbm_media": {"kbase_id": "mbm.media", "kbase_ws": media_workspace},
    "l1_media": {"kbase_id": "l1.media", "kbase_ws": media_workspace},
}

# Load the known growth phenotype data as a pandas DataFrame
growth_phenotype_data = pd.read_csv(
    os.path.join(FILE_DIR, "known_growth_phenotypes.tsv"), sep="\t"
)

# Add a "geneko" ccolumn to the DataFrame and set it to "none"
growth_phenotype_data["geneko"] = "none"

# Add a mediaws column to the DataFrame with the value of "kbase_ws" for
# the media
growth_phenotype_data["mediaws"] = growth_phenotype_data["minimal_media"].apply(lambda x: media[x]["kbase_ws"])

# Add a media column to the DataFrame with the value of "kbase_id" for
# the media
growth_phenotype_data["media"] = growth_phenotype_data["minimal_media"].apply(lambda x: media[x]["kbase_id"])

# Change the title of the "met_id" column to "addtlCpd"
growth_phenotype_data.rename(columns={"met_id": "addtlCpd"}, inplace=True)

# Change the values of the growth column from "Yes", "No", and "Unsure"
# to 1 and 0
growth_phenotype_data["growth"] = growth_phenotype_data["growth"].replace(
    {"Yes": 1, "No": 0, "Unsure": 1}
)

# Save just the columns we need to a new DataFrame
growth_phenotype_data = growth_phenotype_data[
    ["geneko", "mediaws", "media", "addtlCpd", "growth"]
]

# Save the DataFrame to a tab-separated file
growth_phenotype_data.to_csv(
    os.path.join(FILE_DIR, "kbase_phenotype_set.tsv"), sep="\t", index=False
)
