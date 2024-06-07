import os
import unittest
import warnings

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set path to the `test_files` directory
TESTFILE_DIR = os.path.join(os.path.dirname(__file__), "test_files")

# Define a minimal media, but remove any carbon sources. This will be used
# to test that there is no growth when there is no carbon source present. And
# will then be used to test that there is growth when a specific carbon source
# is added to the media (e.g. glucose).
# The reactions listed here must be present in the model, so this may change
# if the model is updated.
# TOOD: Check that this is actually the media composition for what Zac used
ZAC = {
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

# L1 Minimal Media
# A general purpose marine medium for growing coastal algae
# An enriched seawater medium, with everything added to filtered natural seawater
# Does that mean there are other carbon/nitrogen sources in the media?
L1 = {"EX_cpd00007_e0": 20,  # O2_e0
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
      }


class TestGrowthPhenotypes(unittest.TestCase):
    # Check that there is no growth on a media with no carbon sources
    def test_growth_w_0_C(self):
        # Load the model with cobrapy
        model = cobra.io.read_sbml_model("model.xml")

        # Set the media so that there are no carbon sources
        model.medium = ZAC

        # Run the model
        sol = model.optimize()

        # Check that no biomass is produced
        self.assertEqual(sol.objective_value, 0)

    # Test that there is growth on the expected pheontypes, and that there
    # is no growth on the expected phenotypes. Make a plot of the growth
    # phenotypes.
    def test_expected_growth_phenotypes(self):
        # Load the TSV of the growth phenotypes
        growth_phenotypes = pd.read_csv(
            os.path.join(TESTFILE_DIR, "known_growth_phenotypes.tsv"), sep="\t"
        )

        # Load the model
        model = cobra.io.read_sbml_model("model.xml")

        # Loop through the growth phenotpes, and add the carbon source to the
        # minimal media, run FBA and check if the model grows
        pred_growth = []
        for index, row in growth_phenotypes.iterrows():
            # Check if the model has an exchange reaction for the metabolite
            if "EX_" + row["met_id"] + "_e0" in [r.id for r in model.reactions]:
                # If it does, add the exchange reaction to the minimal media used
                minimal_media = eval(row["minimal_media"]).copy()
                minimal_media["EX_" + row["met_id"] + "_e0"] = 1000.0
                # Set the media
                model.medium = minimal_media
                # Run the model
                sol = model.optimize()
                # Check if the model grows
                if sol.objective_value > 0:
                    # If it does, add 'Y' to the list
                    pred_growth.append("Yes")
                    # Give a warning if growth was not expected
                    if row["growth"] == "No":
                        warnings.warn(
                            "Model grew on "
                            + row["c_source"]
                            + ", but growth was not expected."
                        )
                else:
                    # If it doesn't, add 'N' to the list
                    pred_growth.append("No")
                    # Give a warning if growth was expected
                    if row["growth"] == "Yes":
                        warnings.warn(
                            "Model did not grow on "
                            + row["c_source"]
                            + ", but growth was expected."
                        )
            else:
                # If it doesn't have the exchange reaction, add None to the list
                pred_growth.append("No Exchange")
                # Give a warning if growth was expected
                if row["growth"] == "Yes":
                    warnings.warn(
                        "Model did not have an exchange reaction for "
                        + row["c_source"]
                        + ", but growth was expected."
                    )

        # Add the list as a new column in the dataframe
        growth_phenotypes["pred_growth"] = pred_growth

        # Save the dataframe as a TSV
        growth_phenotypes.to_csv(
            os.path.join(TESTFILE_DIR, "known_growth_phenotypes_w_pred.tsv"),
            sep="\t",
            index=False,
        )

        # Plot a categorical heatmap of the growth phenotypes, where the rows
        # are the metabolites and the columns are the experimental and predicted
        # growth phenotypes. Show growth as blue and no growth as orange, and
        # unsure as gray
        # First, make a new dataframe with the metabolites as the rows and the
        # experimental and predicted growth phenotypes as the columns
        # Combine the values of "minimal_media" and "c_source" into one column
        growth_phenotypes["c_source"] = (
            growth_phenotypes["minimal_media"] + " " + growth_phenotypes["c_source"]
        )
        # And set it as the index
        growth_phenotypes = growth_phenotypes.set_index("c_source")
        # Subset the other columns, to have just the growth and predicted growth
        growth_phenotypes = growth_phenotypes[["growth", "pred_growth"]]

        # Rename the columns and the index to be longer/more descriptive
        growth_phenotypes.index.name = "Media/Carbon Source"
        growth_phenotypes = growth_phenotypes.rename(
            columns={
                "growth": "Experimental Growth",
                "pred_growth": "FBA Predicted Growth",
            }
        )

        # Replace all of the "No Exchange" values with "No"
        growth_phenotypes = growth_phenotypes.replace("No Exchange", "No")

        # Make a dictionary for the phenotypes to numbers
        value_to_int = {"Unsure": 0, "No": 1, "Yes": 2}
        n = len(value_to_int)

        # Make a colormap of specified colors (in numerical order for the phenotypes)
        # cmap = ['gray', '#F18F01', '#399E5A'] # Gray, orange, green
        cmap = ["#5E5E5E", "#FF7D0A", "#024064"]  # C-CoMP gray, orange, and dark blue

        # Plot the heatmap
        ax = sns.heatmap(
            growth_phenotypes.replace(value_to_int),
            cmap=cmap,
            linewidths=4,
            linecolor="white",
        )

        # modify colorbar:
        colorbar = ax.collections[0].colorbar
        r = colorbar.vmax - colorbar.vmin
        colorbar.set_ticks([colorbar.vmin + r / n * (0.5 + i) for i in range(n)])
        colorbar.set_ticklabels(list(value_to_int.keys()))

        # Move the x-axis labels to the top
        plt.tick_params(
            axis="both",
            which="major",
            labelsize=10,
            labelbottom=False,
            bottom=False,
            top=False,
            labeltop=True,
        )

        # Make sure that the y-axis labels are not cut off
        plt.tight_layout()

        # Save the figure
        plt.savefig(os.path.join(TESTFILE_DIR, "exp_vs_pred_growth_phenotypes.png"))


if __name__ == "__main__":
    unittest.main()
