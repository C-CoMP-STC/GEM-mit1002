import json
import os
import pickle

import pandas as pd

# This file is for defining the different minimal media with no carbon sources
# for the growth tests. These media are defined as dictionaries, where the keys
# are the exchange reactions for the metabolites in the media, and the values
# are the lower bound for the exchange reaction. The lower bound is set to 1000
# for most of the metabolites, other than oxygen.
# The exchange reactions listed here must be present in the model, so not
# every metabolite in the experimental media can actuually be in the media.
# Where the exchange reaction is missing, I've commented out that line.
# TODO: Instead of setting the media with cobrapy's model.medium(), I could
# write a new function that only sets the exchange reactions that are present
# so that the model doesn't throw an error if the exchange reaction is missing.
# if the model is updated.
# All exchange reactions are defined using the standard modelSEED nomenclature

# minimal_media
# The minimal media I initally used for simualtions/gap filling
# Not necessarily based on anything used in the lab
minimal_media = {
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00007_e0": 20,  # O2_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00075_e0": 1000,  # Nitrite_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00635_e0": 1000,  # Cbl_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# mbm_media
# Minimal Basal Medium used by Zac in the Moran lab for the first round
# of growth tests for A. mac MI1002
# Does not contain an organic carbon source. They added organics to 12 mM
# C-equivalents (e.g. 2 mM glucose = 12 mM C-equivalent)
mbm_media = {
    "EX_cpd00007_e0": 20,  # O2_e0
    # Artificial Sea Water (ASW) Solution
    "EX_cpd00971_e0": 1000,  # Na+_e0 (in NaCl, Na2SO4, NaF, and NaHCO3)
    "EX_cpd00099_e0": 1000,  # Cl-_e0 (in NaCl, MgCl2, CaCl2, and SrCl2)
    "EX_cpd00048_e0": 1000,  # Sulfate (O4S) (in Na2SO4)
    "EX_cpd00205_e0": 1000,  # K+ (in KCl, and KBr)
    "EX_cpd00966_e0": 1000,  # Bromide (Br-) (in KBr)
    "EX_cpd09225_e0": 1000,  # Boric acid (H3BO3) (in H3BO3)
    "EX_cpd00552_e0": 1000,  # Fluoride (F-) (in NaF)
    "EX_cpd00242_e0": 1000,  # Biocarbonate (HCO3-) (in NaHCO3)
    "EX_cpd00254_e0": 1000,  # Mg_e0 (in MgCl2)
    "EX_cpd00063_e0": 1000,  # Ca2+_e0 (in CaCl2)
    "EX_cpd09695_e0": 1000,  # Strontium (Sr2+) (in SrCl2)
    # FeEDTA (Is a chelating agent)
    # TODO: Check if I chould use uncharge Fe instead of Fe3+
    "EX_cpd10516_e0": 1000,  # fe3_e0 (in FeEDTA)
    "EX_cpd00240_e0": 1000,  # EDTA (in FeEDTA)
    # Basal Medium
    "EX_cpd28238_e0": 1000,  # tris-hydrochloride (tris-HCl)
    "EX_cpd00013_e0": 1000,  # Ammonia (in NH4Cl)
    # Cl (from NH4Cl) already included
    # K (from K2HPO4) already included
    "EX_cpd00009_e0": 1000,  # Phosphate (HO4P) (in K2HPO4)
    "EX_cpd00001_e0": 1000,  # H2O (in H2O)
    # Vitamin Supplement
    # Water already included
    "EX_cpd00104_e0": 1000,  # Biotin (Vitamin H)
    "EX_cpd00393_e0": 1000,  # Folate (Folic acid)
    "EX_cpd00263_e0": 1000,  # Pyridoxine (Pyridoxol)
    "EX_cpd00220_e0": 1000,  # Riboflavin
    "EX_cpd00305_e0": 1000,  # Thiamine
    "EX_cpd00133_e0": 1000,  # Nicotinic acid (Niacinamide)
    "EX_cpd00644_e0": 1000,  # Pantothenic acid (Pantothenate)
    "EX_cpd01826_e0": 1000,  # Cyanocobalamin (Dicopac)
    "EX_cpd00443_e0": 1000,  # p-Aminobenzoic acid (ABEE)
    # Not in the media definition, but needed for growth
    "EX_cpd00058_e0": 1000,  # Cu2+_e0  NOT IN MBM MEDIA
    "EX_cpd00030_e0": 1000,  # Mn2+_e0  NOT IN MBM MEDIA
    "EX_cpd00034_e0": 1000,  # Zn2+_e0  NOT IN MBM MEDIA
    "EX_cpd00149_e0": 1000,  # Co2+_e0  NOT IN MBM MEDIA
}

# l1_media
# L1 Minimal Media
# A general purpose marine medium for growing coastal algae
# An enriched seawater medium, with everything added to filtered natural seawater
# FIXME: Does that mean there are other carbon/nitrogen sources in the media?
l1_media = {
    "EX_cpd00007_e0": 20,  # O2_e0
    # L1 salts
    "EX_cpd00971_e0": 1000,  # Na+_e0 (in NaNO3, NaH2PO4, NaSiO3, Na2EDTA, NaMoO4, Na3VO4)
    "EX_cpd00009_e0": 1000,  # Phosphate (HO4P) (in NaH2PO4)
    "EX_cpd20826_e0": 1000,  # Silica (O2Si) (in NaSiO3)
    # Trace element solution
    "EX_cpd00240_e0": 1000,  # EDTA (in Na2EDTA)
    "EX_cpd10516_e0": 1000,  # fe3_e0 (in FeCl3)
    "EX_cpd00099_e0": 1000,  # Cl- (in FeCl3, MnCl2, CoCl2)
    "EX_cpd00030_e0": 1000,  # Mn2+ (in MnCl2)
    "EX_cpd00034_e0": 1000,  # Zn2+ (in ZnSO4)
    "EX_cpd00048_e0": 1000,  # Sulfate (O4S) (in ZnSO4, CuSO4, NiSO4)
    "EX_cpd00149_e0": 1000,  # Co2+ (in CoCl2)
    "EX_cpd00058_e0": 1000,  # Cu2+_e0 (in CuSO4)
    "EX_cpd11574_e0": 1000,  # Molybdate (MoO4) (in NaMoO4)
    "EX_cpd03387_e0": 1000,  # Selenite (O3Se) (in H2SeO3)
    "EX_cpd00244_e0": 1000,  # Ni2+ (in NiSO4)
    "EX_cpd08438_e0": 1000,  # Ortho-vanadate (H2O4V) (in Na3VO4)
    "EX_cpd00205_e0": 1000,  # K+ (in K2CrO4)
    "EX_cpd11595_e0": 1000,  # Chromate (H2CrO4) (in K2CrO4)
    # Vitamin solution
    "EX_cpd00305_e0": 1000,  # Thiamine HCl (Vitamin B1)
    "EX_cpd00104_e0": 1000,  # Biotin (Vitamin H)
    "EX_cpd01826_e0": 1000,  # Cyanocobalamin (Vitamin B12)
    # Not in L1, but needed to grow
    # "EX_cpd00013_e0": 1000,  # NH3_e0 (Needed for growth on glucose and acetate)
    "EX_cpd00063_e0": 1000,  # Ca2+_e0 (Needed for growth on alanine)
    "EX_cpd00254_e0": 1000,  # Mg_e0 (Needed for growth on alanine)
}

# Save a pickle file with the media definitions
with open(os.path.join(os.path.dirname(__file__), "media_definitions.pkl"), "wb") as f:
    pickle.dump(
        {"minimal_media": minimal_media, "mbm_media": mbm_media, "l1_media": l1_media},
        f,
    )

# Need to load in the ModelSEED database first
modelseed_db = json.load(
    open(
        "/Users/helenscott/Documents/PhD/Segre-lab/ModelSEEDDatabase/Biochemistry/compounds.json"
    )
)
# Convert to a dictionary with the ModelSEED IDs as the keys for easier searching
modelseed_db = {met["id"]: met for met in modelseed_db}


# Convert to a TSV file to upload to KBase
def write_media_tsv(media_dict, media_name, modelseed_db):
    media_df = pd.DataFrame.from_dict(media_dict, orient="index", columns=["minFlux"])
    # Fix the names of the compounds
    media_df.index = media_df.index.str.replace("EX_", "").str.replace("_e0", "")
    media_df.index.name = "compounds"
    # Make the min flux negative
    media_df["minFlux"] = -1 * media_df["minFlux"]
    # Set the max flux for everything to be 1000
    media_df["maxFlux"] = 1000
    # Add the names of the compounds
    media_df["name"] = media_df.index.map(lambda x: modelseed_db[x]["name"])
    # Add the formula of the compounds
    media_df["formula"] = media_df.index.map(lambda x: modelseed_db[x]["formula"])
    # Set the concentration to be something?
    # TODO: Does the concentration matter?
    media_df["concentration"] = 1
    # Save
    media_df.to_csv(
        os.path.join(os.path.dirname(__file__), media_name + "_media.tsv"), sep="\t"
    )


write_media_tsv(mbm_media, "mbm", modelseed_db)
write_media_tsv(l1_media, "l1", modelseed_db)


# Write a function to conver the alias strings to a dictionary
def convert_aliases_to_dict(alias_string):
    return {
        alias.split(":")[0]: [ak.strip() for ak in alias.split(":")[1].split(";")]
        for alias in alias_string
        if alias
    }


# Make a new dataframe for the carveme database
media_db = pd.DataFrame(columns=["medium", "description", "compound", "name"])


# Define a function to make the carveme media database
def make_carveme_media(media_dict, media_id, media_name, modelseed_db, media_db):
    for ex_rxn, min_flux in media_dict.items():
        met = ex_rxn.replace("EX_", "").replace("_e0", "")
        name = modelseed_db[met]["name"]
        aliases = convert_aliases_to_dict(modelseed_db[met]["aliases"])
        if "BiGG" not in aliases:
            print(f"No BiGG ID for {name}")
            bigg_to_use = met
        else:
            bigg_id = aliases["BiGG"]
            if len(bigg_id) == 0:
                print(f"No BiGG ID for {name}")
                bigg_to_use = met
            if len(bigg_id) > 1:
                print(f"Multiple BiGG IDs for {name}: {bigg_id}")
            bigg_to_use = bigg_id[0]
        media_db = pd.concat(
            [
                media_db,
                pd.DataFrame(
                    {
                        "medium": media_id,
                        "description": media_name,
                        "compound": bigg_to_use,
                        "name": modelseed_db[met]["name"],
                    },
                    index=[0],
                ),
            ],
            ignore_index=True,
        )
    return media_db


media_db = make_carveme_media(
    mbm_media, "mbm", "Minimal Basal Medium (Moran Lab)", modelseed_db, media_db
)
media_db = make_carveme_media(
    l1_media, "l1", "L1 Minimal Medium", modelseed_db, media_db
)

# Save the media database
media_db.to_csv(
    os.path.join(os.path.dirname(__file__), "no_c_media_database.tsv"),
    sep="\t",
    index=False,
)
