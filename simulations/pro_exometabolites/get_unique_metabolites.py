
import pandas as pd
import os

def get_unique_metabolites():
    """
    Reads the data file and prints the unique metabolite names.
    """
    # Get the directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Data file is in the same directory
    data_file = os.path.join(script_dir, 'ProDiel_quant_20260211.csv')

    # Read the data
    df = pd.read_csv(data_file)

    # Get unique metabolite names
    unique_metabolites = df['CleanName'].unique()

    print("Unique metabolites:")
    for metabolite in unique_metabolites:
        print(metabolite)

if __name__ == '__main__':
    get_unique_metabolites()
