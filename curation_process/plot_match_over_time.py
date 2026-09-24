import os
import sys
import warnings

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd
from adjustText import adjust_text

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(os.path.dirname(__file__))

# Add the repo root to the system path so tools/ is importable
sys.path.append(PROJECT_ROOT)
# Import the shared plot styles
from tools.plot_styles import set_manuscript_style, set_plot_style, summer_colors

# Global figure style (font, sizes, vector text) -- must run before any
# figure or axes is created, see set_manuscript_style's docstring.
set_manuscript_style()

# Define which PRs to highlight on the plot
# These are PRs that caused a significant change in the number of matches or unbounded flux reactions
GROWTH_HIGHLIGHT_PRS = {
    160: "Coneect dead-ends to glyoxylate",
    192: "Removed duplicate thiamine reaction",
    281: "Standardize amino acid transporters (no longer grows on leucine and isoleucine)",
    309: "Add Na+ symporters for leucine and isoleucine",
}
FLUX_HIGHLIGHT_PRS = {
    211: "Remove reactions with glucose anomers",
    344: "Fixed the GAPDHN infeasible loop",
}

# Load the data
data = pd.read_csv(os.path.join(FILE_DIR, "growth_match_summary.csv"))

# Skip rows that say "Error" in the % Match column
data = data[data["% Match"] != "ERROR"]

# Convert columns to numeric type
# TODO: Do I need any other columns numerical? Are there any that can't be?
cols_to_fix = ["Matches", "% Match", "PR Number", "Unbounded Flux Reactions", "Total"]
data[cols_to_fix] = data[cols_to_fix].apply(pd.to_numeric, errors="coerce")

# Sort the data by PR number and reset the index
data = data.sort_values("PR Number").reset_index(drop=True)

# Every row must have been scored the same way, or the line is comparing two
# different definitions of "match". run_tests_on_prs.py re-runs stale rows
# automatically, so a mixture here means the file was edited by hand or a run
# was interrupted part way through.
if "Scoring Version" in data.columns:
    versions = sorted(data["Scoring Version"].dropna().unique())
    if len(versions) > 1:
        raise ValueError(
            f"growth_match_summary.csv mixes scoring versions {versions}; the "
            f"points are not comparable. Re-run run_tests_on_prs.py (it will "
            f"re-evaluate the stale rows)."
        )

# The denominator: conditions with a definite Yes/No observation and no
# exclusion reason. Constant across PRs by construction, which is what makes
# the two ends of the line comparable.
n_interpretable = int(data["Total"].dropna().iloc[-1]) if "Total" in data else None
if n_interpretable is not None and data["Total"].dropna().nunique() > 1:
    warnings.warn(
        "the 'Total' column is not constant across PRs, so the match counts "
        "have different denominators and the line's slope is not meaningful"
    )

# Pick one color per series so the line and its axis can be matched by eye.
# Deliberately not the teal/red of the pipeline schematic in panel A: reviewers
# read the shared colors as a shared meaning, and there isn't one.
matches_color = summer_colors["green"]
flux_color = summer_colors["yellow"]

# Distinct markers as well as distinct colors, so the two series stay separable
# in grayscale and for readers with colour vision deficiency
matches_marker = "o"
flux_marker = "s"

# Two stacked panels sharing one x axis, rather than the twin y axes this used
# to have. The growth-match line gets the taller panel because it is the result
# the figure is about; the unbounded-flux line is a diagnostic, so it gets a
# short panel underneath. Giving each line its own panel also means each y axis
# belongs to exactly one series, so the axes no longer have to be color-coded
# to say which line they describe.
#
# 8.5 in is the full width of the figure as submitted, labels and legends
# included, so nothing gets scaled (and no in-figure text shrinks) afterwards.
fig, (ax_match, ax_flux) = plt.subplots(
    2,
    1,
    sharex=True,
    figsize=(8.5, 5.5),
    height_ratios=[2.5, 1],
)
# Scale the unbounded flux axis to show detail near zero
ax_flux.set_yscale("symlog", linthresh=1)

# Plot the number of matches on the top panel
ax_match.plot(
    data.index,
    data["Matches"],
    marker=matches_marker,
    markersize=4,
    linestyle="-",
    # Slightly heavier than the default: the yellow of this palette is light
    # against white, and a hairline in it disappears at figure scale
    linewidth=1.8,
    color=matches_color,
)

# Plot the number of arbitrarily large reactions on the bottom panel
ax_flux.plot(
    data.index,
    data["Unbounded Flux Reactions"],
    marker=flux_marker,
    # Smaller than the top panel's markers: this panel is a third the height,
    # so markers at the same size merge into a solid slab along the zero line
    markersize=2.5,
    linestyle="-",
    # Slightly heavier than the default: the yellow of this palette is light
    # against white, and a hairline in it disappears at figure scale
    linewidth=1.8,
    color=flux_color,
)

# Apply the shared style (gray axis lines, no top/right spines, gray text)
set_plot_style(ax_match)
set_plot_style(ax_flux)

# Titles and labels
ax_match.set_title("Model Performance Over Time")
if n_interpretable:
    ax_match.set_ylabel(
        f"Growth Phenotypes Matching\nExperimental Data (of {n_interpretable})"
    )
else:
    ax_match.set_ylabel("Growth Phenotypes Matching\nExperimental Data")
# Full range of the metric, so the height of the line reads as a fraction of
# what could be matched rather than being rescaled to whatever was achieved.
# The extra 6% is headroom for the PR labels.
ax_match.set_ylim(0, (n_interpretable or 55) * 1.06)

ax_flux.set_ylabel("Reactions with\nFlux > 100 (Log Scale)")
# Headroom above the tallest spike for that panel's PR labels. Set explicitly
# rather than as a multiple of the data: on a log scale a proportional margin
# is an extra decade or nothing at all depending on where the maximum falls.
# The bottom goes slightly below zero so the stars on the highlighted PRs that
# sit at zero (PR 344) are drawn whole rather than clipped by the x axis. Safe
# on this scale because symlog is linear within +/-linthresh, so a small
# negative bottom is a small amount of space, not a decade.
ax_flux.set_ylim(-0.55, 300)
ax_flux.set_xlabel("Pull Request Number")

# Thin out the x-tick labels: with ~120 points, labeling every PR is unreadable,
# so show every Nth PR number instead. Set on the shared (bottom) axis.
step = max(1, len(data) // 15)
tick_positions = data.index[::step]
ax_flux.set_xticks(tick_positions)
ax_flux.set_xticklabels(data["PR Number"].iloc[::step], rotation=45, ha="right")

# ---- Annotate the highlighted PRs -------------------------------------------
# The x-axis is the row index, not the PR number, so map PR number -> x position
pr_to_index = dict(zip(data["PR Number"], data.index))


def darken(color, factor=0.65):
    """Return a darker shade of a color so the star stands out from its line.

    Same factor for both series. The yellow needs it more than the green does,
    but a per-series factor would make the two stars differ in a way that reads
    as meaningful when it isn't.
    """
    r, g, b = mcolors.to_rgb(color)
    return (r * factor, g * factor, b * factor)


# One left edge for both y-axis labels, instead of each sitting wherever its
# own tick labels end
fig.align_ylabels([ax_match, ax_flux])

# Finalize the layout before annotating so adjustText sees the axes at their
# final size
fig.tight_layout()
fig.canvas.draw()


def add_highlights(ax, y_col, highlights, series_color, expand=(1.4, 1.8)):
    """Star the highlighted PRs on `ax` and lay their labels out around them.

    Each panel is laid out on its own: the labels only ever have to avoid the
    one line drawn in that panel, which is the part that got simpler by
    splitting the twin axes into two panels.
    """
    star_color = darken(series_color)
    label_texts = []
    for pr in highlights:
        if pr not in pr_to_index:
            warnings.warn(f"Highlight PR #{pr} is not in the plotted data; skipping.")
            continue
        x = pr_to_index[pr]
        y = data.loc[x, y_col]
        # Bigger star, darker shade of the line color, white halo to separate it
        ax.plot(
            x,
            y,
            marker="*",
            markersize=22,
            color=star_color,
            markeredgecolor="white",
            markeredgewidth=1.4,
            zorder=7,
        )
        label_texts.append(
            ax.text(
                x,
                y,
                f"PR {pr}",
                color=star_color,
                fontsize=12,
                fontweight="bold",
                ha="center",
                va="center",
                zorder=8,
            )
        )
    # Points for the labels to avoid: every plotted marker of this panel's series
    adjust_text(
        label_texts,
        x=list(data.index),
        y=list(data[y_col]),
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="0.5", lw=0.8),
        expand=expand,
        force_text=(0.4, 0.7),
        ensure_inside_axes=True,
        min_arrow_len=6,
    )


add_highlights(ax_match, "Matches", GROWTH_HIGHLIGHT_PRS, matches_color)
# The flux line sits on zero almost everywhere, so its labels start out on top
# of it and the whole panel above is empty. Push them further vertically than
# the top panel's, where space is tighter and a big shove would strand a label
# far from its star.
add_highlights(
    ax_flux,
    "Unbounded Flux Reactions",
    FLUX_HIGHLIGHT_PRS,
    flux_color,
    expand=(1.4, 4.0),
)

# Save the plot
fig.savefig(os.path.join(FILE_DIR, "match_over_time.png"), dpi=300)
