"""
Compute per-Prochlorococcus-cell metabolite release rates from a diel timecourse.

For each metabolite, computes dC/dt between consecutive timepoints (forward
differences), normalizes by the mean Pro cell density across the interval, and
then  clip negative dC/dt to zero. This treats intervals
where Pro is reabsorbing a metabolite as supplying zero flux to a
heterotrophic neighbor.

The output rate at row i is the rate for the interval [t_i, t_{i+1}], assigned
to the start of the interval. Units are nmol per Pro cell per hour. This is
NOT yet normalized to Alteromonas dry weight; that conversion happens later
when these rates are turned into FBA exchange bounds.

Inputs:
    ProDiel_filtered_meanByTimepoint.csv: filtered exometabolite concentrations
        (long format: CleanName, timepoint, mean_nM, sd_nM, n_reps)
    extrapolated_cellcounts.csv: Pro cell density timecourse from Natalie

Output:
    ProDiel_per_pro_cell_rates.csv: long format with columns
        metabolite, interval_start_h, interval_end_h, dC_dt_nM_per_hr,
        mean_pro_density_cells_per_mL, per_pro_cell_rate_nmol_per_hr,
        per_pro_cell_rate_nmol_per_hr_raw, cell_density_quality
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Constants for file paths
FILE_DIR = Path(__file__).parent
OUT_DIR = Path(FILE_DIR, "results")
# If the output directory doesn't exist, create it
OUT_DIR.mkdir(exist_ok=True)

# Set the paths/names for the input and output files
INPUT_METABOLITES = Path(FILE_DIR, "data", "ProDiel_filtered_meanByTimepoint.csv")
INPUT_CELLCOUNTS = Path(FILE_DIR, "data", "extrapolated_cellcounts.csv")
OUTPUT_FILE = Path(OUT_DIR, "ProDiel_per_pro_cell_rates.csv")
PLOT_FILE = Path(OUT_DIR, "ProDiel_per_pro_cell_rates.png")

# Define the experiments dark periods (from Katie's chapter Figure 4.3 shading)
DARK_PERIODS = [(10, 22), (34, 46)]


def main() -> None:
    # Load the data
    metabolites = load_metabolites(INPUT_METABOLITES)
    cell_densities = load_cell_densities(INPUT_CELLCOUNTS)

    # Calculate per-Pro-cell rates
    rates = compute_per_pro_cell_rates(metabolites, cell_densities)

    # Save the results
    rates.to_csv(OUTPUT_FILE, index=False)

    # Plot the results (optional, for sanity checking)
    plot_rates(rates, PLOT_FILE)


def load_cell_densities(path: Path) -> pd.DataFrame:
    """Load cell density data and deduplicate to one row per timepoint.

    cell_count_mean is the same across the three replicates at each timepoint,
    so we keep just one representative row per time.
    """
    df = pd.read_csv(path)
    dedup = (
        df.drop_duplicates("time_h")[
            ["time_h", "cell_count_mean", "cell_count_sd", "typeOfData"]
        ]
        .sort_values("time_h")
        .reset_index(drop=True)
    )
    return dedup


def load_metabolites(path: Path) -> pd.DataFrame:
    """Load filtered metabolite concentrations."""
    df = pd.read_csv(path)
    return df.sort_values(["CleanName", "timepoint"]).reset_index(drop=True)


def compute_per_pro_cell_rates(
    metabolites: pd.DataFrame,
    cell_densities: pd.DataFrame,
) -> pd.DataFrame:
    """Compute per-Pro-cell release rates via forward differences.

    For each consecutive pair of timepoints (t_i, t_{i+1}) at which a metabolite
    is measured, computes:
        dC/dt        = (C(t_{i+1}) - C(t_i)) / (t_{i+1} - t_i)   [nM/hr]
        N_mean       = (N(t_i) + N(t_{i+1})) / 2                  [cells/mL]
        rate_raw     = dC_dt / (N_mean * 1000)                    [nmol/cell/hr]
                       (the 1000 converts cells/mL -> cells/L)
        rate_clipped = max(rate_raw, 0)                            [Option A]
    """
    rows = []
    cd_lookup = cell_densities.set_index("time_h")

    for metabolite, group in metabolites.groupby("CleanName"):
        g = group.sort_values("timepoint").reset_index(drop=True)

        for i in range(len(g) - 1):
            t_start = g.loc[i, "timepoint"]
            t_end = g.loc[i + 1, "timepoint"]
            dt = t_end - t_start

            c_start_nM = g.loc[i, "mean_nM"]
            c_end_nM = g.loc[i + 1, "mean_nM"]
            dC_dt_nM_per_hr = (c_end_nM - c_start_nM) / dt

            # Look up cell densities at both endpoints
            if t_start not in cd_lookup.index or t_end not in cd_lookup.index:
                continue

            n_start = cd_lookup.loc[t_start, "cell_count_mean"]
            n_end = cd_lookup.loc[t_end, "cell_count_mean"]
            type_start = cd_lookup.loc[t_start, "typeOfData"]
            type_end = cd_lookup.loc[t_end, "typeOfData"]

            mean_density_cells_per_mL = (n_start + n_end) / 2
            mean_density_cells_per_L = mean_density_cells_per_mL * 1000

            # nM/hr = nmol/L/hr; divide by cells/L -> nmol/(cell*hr)
            rate_raw = dC_dt_nM_per_hr / mean_density_cells_per_L
            rate_clipped = max(rate_raw, 0)

            quality = (
                "measured"
                if type_start == "measured" and type_end == "measured"
                else "partial_estimate"
            )

            rows.append(
                {
                    "metabolite": metabolite,
                    "interval_start_h": t_start,
                    "interval_end_h": t_end,
                    "dC_dt_nM_per_hr": dC_dt_nM_per_hr,
                    "mean_pro_density_cells_per_mL": mean_density_cells_per_mL,
                    "per_pro_cell_rate_nmol_per_hr": rate_clipped,
                    "per_pro_cell_rate_nmol_per_hr_raw": rate_raw,
                    "cell_density_quality": quality,
                }
            )

    return pd.DataFrame(rows)


def plot_rates(rates: pd.DataFrame, output_path: Path) -> None:
    metabolites = sorted(rates["metabolite"].unique())
    n = len(metabolites)
    ncols = 3
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 2.4 * nrows), sharex=True)
    axes = axes.flatten()

    for ax, met in zip(axes, metabolites):
        sub = rates[rates["metabolite"] == met].sort_values("interval_start_h")
        # Plot at the midpoint of each interval; width matches the interval
        midpoints = (sub["interval_start_h"] + sub["interval_end_h"]) / 2
        widths = sub["interval_end_h"] - sub["interval_start_h"]
        values = sub["per_pro_cell_rate_nmol_per_hr_raw"]

        colors = ["#2c7fb8" if v > 0 else "#d95f0e" for v in values]
        ax.bar(
            midpoints,
            values,
            width=widths * 0.9,
            color=colors,
            edgecolor="black",
            linewidth=0.3,
        )

        # Mark intervals with partial-estimate cell density (hatching)
        partial = sub["cell_density_quality"] == "partial_estimate"
        if partial.any():
            ax.bar(
                midpoints[partial],
                values[partial],
                width=widths[partial] * 0.9,
                facecolor="none",
                edgecolor="black",
                hatch="///",
                linewidth=0.3,
            )

        # Diel shading
        ymin, ymax = ax.get_ylim()
        for d_start, d_end in DARK_PERIODS:
            ax.axvspan(d_start, d_end, color="gray", alpha=0.15, zorder=0)

        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_title(met, fontsize=9)
        ax.tick_params(labelsize=8)
        ax.set_xlim(-1, 47)

    # Hide any unused axes
    for ax in axes[n:]:
        ax.set_visible(False)

    # Shared labels
    fig.supxlabel("Time (h)", fontsize=11)
    fig.supylabel("Per-Pro-cell rate (nmol/cell/hr)", fontsize=11)

    # Legend
    from matplotlib.patches import Patch

    legend_handles = [
        Patch(
            facecolor="#2c7fb8", edgecolor="black", label="Net release (used by FBA)"
        ),
        Patch(
            facecolor="#d95f0e",
            edgecolor="black",
            label="Net reabsorption (clipped to 0)",
        ),
        Patch(
            facecolor="white",
            edgecolor="black",
            hatch="///",
            label="Partial-estimate cell density",
        ),
        Patch(facecolor="gray", alpha=0.3, label="Dark period"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=4,
        fontsize=9,
        bbox_to_anchor=(0.5, -0.01),
    )

    plt.tight_layout(rect=(0.02, 0.04, 1, 1))
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
