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

# Load the excel sheet with the metabolite groups and order I want to use
metabolite_groups = pd.read_excel(input_path, sheet_name="4.Nomenclature_mapping", skiprows=1)

# Split the table into dataframes for each group by taking the rows between the group names
groups = ["PROTEIN", "tRNA", "mRNA and DNA", "LIPIDS", "COFACTORS & PROSTHETIC GROUPS", "CARBOHYDRATES, CELL-WALL and OTHERS"]
dfs = []
for i, group in enumerate(groups):
    # Get the index of where that group will start in the table (including the row with the group name)
    start_row = metabolite_groups.index[metabolite_groups['This work'] == group].tolist()[0]
    if i < len(groups) - 1:
        # Get the index of where the next group will start in the table (including the row with the group name)
        end_row = metabolite_groups.index[metabolite_groups['This work'] == groups[i + 1]].tolist()[0]
        # Get the rows between the group names
        df_group = metabolite_groups.loc[start_row+1:end_row-1]
    else:
        # For the last group, the group will end at the end of the table
        end_row = metabolite_groups.shape[0]
        # Get the rows between the group name and the end of the table
        df_group = metabolite_groups.loc[start_row+1:end_row]
    # Add a column with the group name
    df_group['Group'] = group
    # Get rid of all the columns other than "This work" and "Group", and rename "This work" to "Metabolite"
    df_group = df_group[['This work', 'Group']].rename(columns={'This work': 'Metabolite'})
    # Add to the list of dataframes
    dfs.append(df_group)

# Merge all of the dataframes into a single dataframe
group_metabolite_ids = pd.concat(dfs)

# Add a row called 'Group' to the pivot table, which will contain the metabolite group for each metabolite
df_pivot['Group'] = df_pivot.index.map(group_metabolite_ids.set_index('Metabolite')['Group'])

# Reorder to match the order of the metabolite groups
df_pivot = df_pivot.sort_values('Group')

# Save the result to a CSV file
output_path = os.path.join(OUT_DIR, 'standardized_compositions.csv')
df_pivot.to_csv(output_path)
