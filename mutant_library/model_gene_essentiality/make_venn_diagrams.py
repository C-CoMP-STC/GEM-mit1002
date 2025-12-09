import os

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2, venn3

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(FILE_DIR, "results")
MUTANT_LIB_DIR = os.path.dirname(FILE_DIR)

########################################################################
# Venn Diagram of the gene locus tags from the genomes from the mutant
# library and the model
########################################################################
# Load the data from the mutant library lookup file
lookup_file_path = os.path.join(
    MUTANT_LIB_DIR, "MIT1002_mutant-library-gene-lookup.csv"
)
df = pd.read_csv(lookup_file_path)

# Define the conditions for each set
# A gene is considered present if its corresponding cell is not empty (not NaN)
in_mut_lib = df["locus_tag_from_mut_lib"].notna()
in_model = df["locus_tag_from_model"].notna()

# Calculate the sizes of the different sections of the Venn diagram
# 1. Intersection: Present in both columns
intersection_count = (in_mut_lib & in_model).sum()

# 2. Only in Mutant Library: Present in the first column but NOT the second
mut_lib_only_count = (in_mut_lib & ~in_model).sum()

# 3. Only in Model: Present in the second column but NOT the first
model_only_count = (in_model & ~in_mut_lib).sum()

# Create the Venn diagram
plt.figure(figsize=(10, 7))  # Set the figure size for better readability

# venn2 function takes a tuple of (Ab, aB, AB)
# Ab: Size of set A not in B
# aB: Size of set B not in A
# AB: Size of intersection
v = venn2(
    subsets=(mut_lib_only_count, model_only_count, intersection_count),
    set_labels=("Mutant Library", "Model"),
)

# Customize the plot
plt.title("Venn Diagram of Gene Locus Tags", fontsize=16)

# Save the figure to a file
plt.savefig(
    os.path.join(OUT_DIR, "genes_in_model_vs_mutant_lib.png"), bbox_inches="tight"
)

########################################################################
# Venn Diagram of the genes essential in the model and present in the
# mutant library
########################################################################
# Load the data from the merged essentiality file
essentiality_file_path = os.path.join(OUT_DIR, "merged_ko_results_with_mut_lib_id.csv")
df_essentiality = pd.read_csv(essentiality_file_path)

# Define the three sets based on conditions
is_essential = df_essentiality["essential_in_model"] == True
is_not_essential = df_essentiality["essential_in_model"] == False
in_mutant_library = df_essentiality["barcode"].notna()

# Calculate the sizes of the 7 sections for a three-circle Venn diagram
# Let A = Essential, B = Not Essential, C = In Mutant Library

# A and B are mutually exclusive, so any intersection between them is 0.
# A and B, but not C
essential_and_not_essential_only = 0
# A and B and C
all_three = 0

# A only (Essential, but not in Mutant Library)
essential_only = (is_essential & ~in_mutant_library).sum()

# B only (Not Essential, but not in Mutant Library)
not_essential_only = (is_not_essential & ~in_mutant_library).sum()

# C only (In Mutant Library, but has no growth data in the model)
# This can happen if a gene from the library isn't in the model at all.
in_mutant_library_only = (in_mutant_library & ~is_essential & ~is_not_essential).sum()

# Intersection of A and C (Essential AND in Mutant Library)
essential_and_in_mutant_library = (is_essential & in_mutant_library).sum()

# Intersection of B and C (Not Essential AND in Mutant Library)
not_essential_and_in_mutant_library = (is_not_essential & in_mutant_library).sum()


# Create the three-circle Venn diagram
plt.figure(figsize=(12, 8))

# The venn3 function takes a 7-element tuple in a specific order:
# (Abc, aBc, ABc, abC, AbC, aBC, ABC)
# Our mapping: A=Essential, B=Not Essential, C=In Mutant Library
# Since A & B are mutually exclusive, ABc and ABC are 0.
# We use a dummy placeholder for the B set in the function call.
venn3(
    subsets=(
        essential_only,  # A, not C
        in_mutant_library_only,  # C, not A or B
        essential_and_in_mutant_library,  # A and C
        not_essential_only,  # B, not C
        0,  # B and A (impossible)
        not_essential_and_in_mutant_library,  # B and C
        0,  # A and B and C (impossible)
    ),
    set_labels=("Essential in Model", "In Mutant Library", "Not Essential in Model"),
)

plt.title(
    "Comparison of Model Gene Essentiality and Mutant Library Coverage", fontsize=16
)

# Save the figure
plt.savefig(
    os.path.join(OUT_DIR, "essentiality_vs_mutant_lib_venn.png"), bbox_inches="tight"
)
