import os
import pandas as pd

# Set the output directory (where the results.pkl file will be saved)
OUT_DIR = os.path.dirname(os.path.realpath(__file__))

# Load the Excel file containing the composition table
input_path = os.path.join(OUT_DIR, 'xavier_2017_si.xlsx')
# Read the table from sheet number 2, but skip the first row (header)
# Requires the openpyxl package to be installed
# Works on the "gem-reconstruction" conda environment
df = pd.read_excel(input_path, sheet_name="2.Standardized_nom_compositions", skiprows=1)

# Clean the DataFrame by removing the unnecessary metadata rows and keeping the relevant data
# We start extracting from row 3 (index 2) and rename the columns appropriately
df_clean = df.iloc[2:, :].reset_index(drop=True)

# Set the column names to the organism names given in row 1 (index 1)
df_clean.columns = df.iloc[1, :].values

# Remove the first column (which is currently empty)
df_clean = df_clean.drop(columns=df_clean.columns[0])

# Melt the DataFrame to have metabolites in a single column and organism names as another column
df_melted = df_clean.melt(var_name='Model', value_name='Metabolite')

# Drop rows with NaN values in the 'Metabolite' column
df_melted.dropna(subset=['Metabolite'], inplace=True)

# Create a pivot table where each column is a model and each row is a metabolite, indicating presence (1) or absence (0)
df_pivot = pd.pivot_table(df_melted, index='Metabolite', columns='Model', aggfunc='size', fill_value=0)

# Save the result to a CSV file
output_path = os.path.join(OUT_DIR, 'standardized_compositions.csv')
df_pivot.to_csv(output_path)
