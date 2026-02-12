import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_data(df, in_or_ex, output_dir):
    """
    Plots the data for a given category (IN or EX).
    """
    g = sns.relplot(
        data=df,
        x="timepoint",
        y="nM",
        hue="CleanName",
        style="replicateID",
        kind="line",
        col="CleanName",
        col_wrap=4,
        facet_kws=dict(sharey=False),
    )
    g.fig.suptitle(f"{in_or_ex} Metabolites", y=1.03)
    plt.savefig(os.path.join(output_dir, f"{in_or_ex}_metabolites.png"))
    plt.close()


def main():
    """
    Main function to read data and generate plots.
    """
    # Get the directory of the script
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Data file is in the same directory
    data_file = os.path.join(script_dir, "ProDiel_quant_20260211.csv")

    # Create an output directory for plots
    output_dir = os.path.join(script_dir, "plots")
    os.makedirs(output_dir, exist_ok=True)

    # Read the data
    df = pd.read_csv(data_file)

    # Separate data into IN and EX
    df_ex = df[df["INorEX"] == "EX"]
    df_in = df[df["INorEX"] == "IN"]

    # Plot EX data
    if not df_ex.empty:
        print("Plotting EX data...")
        plot_data(df_ex, "EX", output_dir)
        print(f"Saved EX plot to {os.path.join(output_dir, 'EX_metabolites.png')}")

    # Plot IN data
    if not df_in.empty:
        print("Plotting IN data...")
        plot_data(df_in, "IN", output_dir)
        print(f"Saved IN plot to {os.path.join(output_dir, 'IN_metabolites.png')}")


if __name__ == "__main__":
    main()
