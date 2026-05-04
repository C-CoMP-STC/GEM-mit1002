"""
Mapping from iSO595v6 (Prochlorococcus) metabolite IDs to ModelSEED IDs
used by the Alteromonas MIT1002 model (model.xml).

Apply this mapping when setting up co-culture COMETS simulations so that
shared extracellular metabolites are recognised as the same compound.
The published iSO595v6.xml is never modified — renaming happens in memory.

Usage
-----
    import cobra
    import cometspy as c
    from pro_met_id_mapping import rename_pro_metabolites

    pro_cobra = cobra.io.read_sbml_model("path/to/iSO595v6.xml")
    rename_pro_metabolites(pro_cobra)
    pro = c.model(pro_cobra)

Notes
-----
- Metabolites present in Pro but absent from Alt (no mapping available):
    HCO3[e]            bicarbonate — not in Alt model
    Glycolate[e]       key Pro exudate — not yet in Alt model (consider adding)
    S_Malate[e]        malate — not in Alt model
    Molybdenum[e]      molybdate — not in Alt model
    Selenate[e]        selenate — not in Alt model
    Strontium_cation[e] strontium — not in Alt model
    Cadmium[e]         cadmium — not in Alt model
    Photon[e]          light (Pro-specific, not applicable to Alt)
"""

# Maps Pro cobra metabolite ID -> Alt ModelSEED metabolite ID
PRO_TO_MODELSEED: dict[str, str] = {
    # Inorganic / ions
    "Ammonia[e]": "cpd00013_e0",  # NH3
    "CO2[e]": "cpd00011_e0",  # CO2
    "H2O[e]": "cpd00001_e0",  # H2O
    "H[e]": "cpd00067_e0",  # H+
    "Oxygen[e]": "cpd00007_e0",  # O2
    "Orthophosphate[e]": "cpd00009_e0",  # Phosphate
    "Sulfate[e]": "cpd00048_e0",  # Sulfate
    "Magnesium_cation[e]": "cpd00254_e0",  # Mg2+
    "K[e]": "cpd00205_e0",  # K+
    "Sodium_cation[e]": "cpd00971_e0",  # Na+
    "Calcium_cation[e]": "cpd00063_e0",  # Ca2+
    "Chloride_ion[e]": "cpd00099_e0",  # Cl-
    "Cobalt_ion[e]": "cpd00149_e0",  # Co2+
    "Copper[e]": "cpd00058_e0",  # Cu2+
    "Zn2[e]": "cpd00034_e0",  # Zn2+
    "Fe2[e]": "cpd10516_e0",  # Fe3+ (closest in Alt)
    "Hydrogen_sulfide[e]": "cpd00239_e0",  # HS- / bisulfide
    "Thiosulfate[e]": "cpd00268_e0",  # S2O3
    # Carbon compounds
    "Acetate[e]": "cpd00029_e0",
    "D_Glucose[e]": "cpd00027_e0",
    "Pyruvate[e]": "cpd00020_e0",
    "Succinate[e]": "cpd00036_e0",
    "Fumarate[e]": "cpd00106_e0",
    "Ethanol[e]": "cpd00363_e0",
    "Formate[e]": "cpd00047_e0",
    "CO[e]": "cpd00204_e0",
    "Methanol[e]": "cpd00116_e0",
    "Citrate[e]": None,  # not in Alt
    # Amino acids
    "Glycine[e]": "cpd00033_e0",
    "L_Alanine[e]": "cpd00035_e0",
    "L_Arginine[e]": "cpd00051_e0",
    "L_Asparagine[e]": "cpd00132_e0",
    "L_Aspartate[e]": None,  # not in Alt
    "L_Glutamate[e]": "cpd00023_e0",
    "L_Glutamine[e]": None,  # not in Alt
    "L_Histidine[e]": None,  # not in Alt
    "L_Isoleucine[e]": "cpd00322_e0",
    "L_Leucine[e]": "cpd00107_e0",
    "L_Lysine[e]": "cpd00039_e0",
    "L_Methionine[e]": None,  # not in Alt
    "L_Phenylalanine[e]": None,  # not in Alt
    "L_Proline[e]": "cpd00129_e0",
    "L_Serine[e]": "cpd00054_e0",
    "L_Threonine[e]": "cpd00161_e0",
    "L_Tryptophan[e]": None,  # not in Alt
    "L_Tyrosine[e]": "cpd00069_e0",
    "L_Valine[e]": "cpd00156_e0",
    # Nucleobases / nucleosides / vitamins
    "Adenine[e]": "cpd00128_e0",
    "Guanosine[e]": "cpd00311_e0",
    "Thymidine[e]": "cpd00184_e0",
    "Cobamide_coenzyme[e]": "cpd03424_e0",  # Vitamin B12
}

# Inverse mapping for convenience
MODELSEED_TO_PRO: dict[str, str] = {
    v: k for k, v in PRO_TO_MODELSEED.items() if v is not None
}


def rename_pro_metabolites(pro_cobra_model) -> None:
    """Rename Pro extracellular metabolite IDs to ModelSEED IDs in-place.

    Modifies the cobra model object — does NOT write to the XML file.
    Metabolites without a ModelSEED mapping are left unchanged.

    Parameters
    ----------
    pro_cobra_model : cobra.Model
        A cobra model loaded from iSO595v6.xml.
    """
    renamed, skipped = [], []
    for met in pro_cobra_model.metabolites:
        new_id = PRO_TO_MODELSEED.get(met.id)
        if new_id is not None:
            met.id = new_id
            renamed.append(new_id)
        else:
            skipped.append(met.id)

    print(
        f"rename_pro_metabolites: renamed {len(renamed)}, left unchanged {len(skipped)}"
    )
