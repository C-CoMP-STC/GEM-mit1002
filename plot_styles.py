import matplotlib.pyplot as plt

# Define the colors from the C-CoMP pallette
LIGHT_BLUE = "#3CB3C0"
DARK_BLUE = "#024064"
ORANGE = "#FF6C2C"
# Define light orange, used on ISME poster
LIGHT_ORANGE = "#FFBB62"


# Define the style for the plots (gray axes, no top or right axis lines)
def set_plot_style(g):
    # Make the axis lines gray
    g.spines["bottom"].set_color("gray")
    g.spines["left"].set_color("gray")
    # Make the tick marks gray
    g.tick_params(axis="x", colors="gray")
    g.tick_params(axis="y", colors="gray")
    # Remove the top and right axis lines
    g.spines["top"].set_visible(False)
    g.spines["right"].set_visible(False)


def carbon_fates_bar(data):
    # Check that the column names are correct
    assert set(data.columns) == set(["co2", "organic_c", "biomass"])
    # Set the column order
    data = data[["biomass", "organic_c", "co2"]]
    # Plot the stacked bar plot
    g = data.plot(kind="bar", stacked=True, color=[DARK_BLUE, LIGHT_BLUE, LIGHT_ORANGE])
    # Move the legend outside of the plot
    custom_labels = ["Biomass", "Organic C", "CO2"]
    lgd = plt.legend(
        bbox_to_anchor=(1.25, 0.5),
        loc="center right",
        borderaxespad=0.0,
        ncol=1,
        labels=custom_labels,
    )
    # Make the legend text gray too
    for text in lgd.get_texts():
        text.set_color("gray")
    # Adjust the bottom margin
    plt.subplots_adjust(bottom=0.2)
    # Style
    set_plot_style(g)
    # Title the plot and make it gray
    g.set_title("Fate of Carbon", color="gray")

    # Return the plot
    return g
