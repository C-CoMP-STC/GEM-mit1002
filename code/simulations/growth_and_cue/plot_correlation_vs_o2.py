"""How the growth rate-CUE correlation across substrates changes with O2.

`plot_cue_vs_growth.py` shows the substrate panel in the growth-rate/CUE
plane: a dot per substrate at the highest O2 level (where the correlation is
essentially perfect and positive) and arrows tracking each substrate as O2
falls. This script collapses that picture to one number per O2 level -- the
Pearson correlation between growth rate and CUE *across* substrates -- and
plots it against O2, so the direction and strength of the relationship can be
read off directly instead of inferred from the arrow fan.

Two figures are produced, one per style of O2 limitation run by
run_growth_cue.py:

  bounds       every substrate gets the same O2 uptake lower bound
               (O2_SET_LEVELS), so at a given level the substrates are not
               equally limited -- one that needs little O2 may still be
               saturated while another is starved.
  percentiles  every substrate gets the same *fraction* of its own saturating
               O2 uptake (O2_PERCENTAGE_LEVELS), so "50%" means the same
               degree of limitation for all of them.

Spread is a percentile bootstrap over substrates: the 22 substrates are
resampled with replacement, r is recomputed, and the 2.5th/97.5th percentiles
are drawn as a band. It answers "how much does this r depend on which
substrates happen to be in the panel", which is the relevant uncertainty here
(the FBA solutions themselves are deterministic, so there is no within-point
error to show). The band is asymmetric near |r| = 1, as it should be.

Outputs
-------
figures/growth_cue_correlation_vs_o2_bounds.{png,svg}
figures/growth_cue_correlation_vs_o2_percentiles.{png,svg}
results/growth_cue_correlation_vs_o2.csv   (r, p, CI and n for every level)
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

FILE_PATH = Path(__file__).resolve().parent
IN_PATH = FILE_PATH / "results"
OUT_PATH = FILE_PATH / "figures"
OUT_PATH.mkdir(exist_ok=True)

sys.path.insert(0, str(FILE_PATH.parents[1]))  # make `tools` importable
from tools.plot_styles import set_manuscript_style, set_plot_style, summer_colors

# Global figure style (font, sizes, vector text) -- must run before any figure
# or axes is created, see set_manuscript_style's docstring.
set_manuscript_style()

# The O2 levels to report, per style. These must stay in step with
# O2_SET_LEVELS and O2_PERCENTAGE_LEVELS in run_growth_cue.py; rather than
# importing them from there (which would pull in cobra and re-run nothing
# useful), the levels are listed here and checked against the results file by
# select_level(), which raises if a level is missing. Levels that produced no
# growth for any substrate (0% O2) are never written to the results file, so
# they are dropped with a warning instead of raising.
BOUND_LEVELS = [1, 10, 20, 30, 40, 50]
PERCENT_LEVELS = list(range(0, 110, 10))

# Substrates needed at a level before a correlation is reported. Low O2 levels
# drop substrates that stop growing, and r over a handful of points is not
# worth plotting.
MIN_SUBSTRATES = 6

# Bootstrap settings. The seed is fixed so re-running gives byte-identical
# figures and CSV, the same reason FLUX_DECIMALS exists in run_growth_cue.py.
N_BOOTSTRAP = 10_000
BOOTSTRAP_SEED = 0
CI_PERCENTILES = (2.5, 97.5)

# Float tolerance for matching a requested level against the values in the
# results file (o2_percent is computed as level / saturation * 100, so "30" can
# land as 30.000000000000004).
LEVEL_TOL = 1e-6

LINE_C = summer_colors["teal"]
BAND_C = summer_colors["teal"]
ZERO_C = "#999999"


def main():
    summary_df = pd.read_csv(IN_PATH / "growth_and_cue.csv")

    bounds_df = correlation_table(summary_df, "o2_bound", BOUND_LEVELS)
    percent_df = correlation_table(summary_df, "o2_percent", PERCENT_LEVELS)

    plot_correlation_vs_o2(
        bounds_df,
        OUT_PATH,
        "growth_cue_correlation_vs_o2_bounds",
        xlabel="O$_2$ uptake bound (mmol O$_2$ gDW$^{-1}$ h$^{-1}$)",
        subtitle="Same O$_2$ bound for every substrate",
    )
    plot_correlation_vs_o2(
        percent_df,
        OUT_PATH,
        "growth_cue_correlation_vs_o2_percentiles",
        xlabel="O$_2$ supply (% of each substrate's saturating uptake)",
        subtitle="Same fraction of saturating O$_2$ for every substrate",
    )

    # One tidy table for both styles, so the numbers in the figures are
    # traceable without re-running the bootstrap.
    out_df = pd.concat(
        [bounds_df.assign(o2_style="bounds"), percent_df.assign(o2_style="percent")],
        ignore_index=True,
    )
    out_df = out_df[
        ["o2_style", "o2_level", "n_substrates", "r", "p", "ci_low", "ci_high"]
    ]
    for col in ["r", "ci_low", "ci_high"]:
        out_df[col] = out_df[col].round(4)
    out_df.to_csv(IN_PATH / "growth_cue_correlation_vs_o2.csv", index=False)
    print(out_df.to_string(index=False))


def select_level(df: pd.DataFrame, level_col: str, level: float) -> pd.DataFrame:
    """Rows at one O2 level, matched on `level_col` with a float tolerance."""
    return df[np.isclose(df[level_col], level, atol=LEVEL_TOL, rtol=0)]


def correlation_table(
    summary_df: pd.DataFrame, level_col: str, levels: list
) -> pd.DataFrame:
    """Pearson r (with bootstrap CI) of growth rate vs CUE at each O2 level.

    One row per level in `levels` that the results file actually contains and
    that has at least MIN_SUBSTRATES growing substrates."""
    rows = []
    for level in levels:
        sub = select_level(summary_df, level_col, level)
        # Only substrates that grew contribute; CUE is undefined without uptake
        sub = sub.dropna(subset=["growth_rate", "cue"])
        if sub.empty:
            print(f"  WARN  no results at {level_col} = {level:g}; skipping")
            continue
        if sub["substrate"].duplicated().any():
            dupes = sorted(sub.loc[sub["substrate"].duplicated(), "substrate"])
            raise ValueError(
                f"Multiple rows for the same substrate at {level_col} = {level:g}: "
                f"{dupes}. The O2 levels in this script no longer match the ones "
                "run by run_growth_cue.py."
            )
        if len(sub) < MIN_SUBSTRATES:
            print(
                f"  WARN  only {len(sub)} substrates at {level_col} = {level:g} "
                f"(< {MIN_SUBSTRATES}); skipping"
            )
            continue

        x = sub["growth_rate"].to_numpy(float)
        y = sub["cue"].to_numpy(float)
        r, p = stats.pearsonr(x, y)
        ci_low, ci_high = bootstrap_ci(x, y)
        rows.append(
            {
                "o2_level": float(level),
                "n_substrates": len(sub),
                "r": r,
                "p": p,
                "ci_low": ci_low,
                "ci_high": ci_high,
            }
        )

    if not rows:
        raise ValueError(f"No usable O2 levels found in the '{level_col}' column")
    return pd.DataFrame(rows).sort_values("o2_level").reset_index(drop=True)


def bootstrap_ci(x, y, n_boot=N_BOOTSTRAP, seed=BOOTSTRAP_SEED):
    """Percentile bootstrap CI for Pearson r, resampling substrates.

    Resamples that happen to be constant in x or y have an undefined r and are
    dropped rather than counted as zero."""
    rng = np.random.default_rng(seed)
    n = len(x)
    idx = rng.integers(0, n, size=(n_boot, n))
    xs, ys = x[idx], y[idx]
    xs = xs - xs.mean(axis=1, keepdims=True)
    ys = ys - ys.mean(axis=1, keepdims=True)
    denom = np.sqrt((xs**2).sum(axis=1) * (ys**2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        rs = (xs * ys).sum(axis=1) / denom
    rs = rs[np.isfinite(rs)]
    if rs.size == 0:
        return np.nan, np.nan
    return tuple(np.percentile(rs, CI_PERCENTILES))


def plot_correlation_vs_o2(
    corr_df: pd.DataFrame,
    out_path: Path,
    filename: str,
    xlabel: str,
    subtitle: str,
    alpha=0.05,
) -> None:
    """Line plot of r vs O2 level with the bootstrap band.

    Filled markers are significant at `alpha`, open markers are not, so the
    strength of the correlation and the confidence in its sign are both
    readable without a second panel."""
    x = corr_df["o2_level"].to_numpy(float)
    r = corr_df["r"].to_numpy(float)

    fig, ax = plt.subplots(figsize=(5.5, 4))

    # r = 0 reference, so the sign of the correlation is readable at a glance
    ax.axhline(0, color=ZERO_C, lw=0.8, ls="--", zorder=1)

    ax.fill_between(
        x,
        corr_df["ci_low"],
        corr_df["ci_high"],
        color=BAND_C,
        alpha=0.18,
        linewidth=0,
        zorder=2,
    )
    ax.plot(x, r, color=LINE_C, lw=2, zorder=3)

    # Significant points get a filled marker, non-significant ones a hollow one
    sig = corr_df["p"].to_numpy(float) < alpha
    ax.scatter(
        x[sig], r[sig], s=45, color=LINE_C, edgecolor="white", linewidth=0.8, zorder=4
    )
    ax.scatter(
        x[~sig], r[~sig], s=45, facecolor="white", edgecolor=LINE_C, lw=1.5, zorder=4
    )

    # Flag levels where substrates dropped out (no growth), since r there is
    # over a different, smaller panel than the rest of the line
    n_max = int(corr_df["n_substrates"].max())
    for xi, ri, ni in zip(x, r, corr_df["n_substrates"]):
        if ni < n_max:
            ax.annotate(
                f"n = {ni}",
                (xi, ri),
                xytext=(0, -12),
                textcoords="offset points",
                fontsize=7,
                ha="center",
                va="top",
                color="gray",
            )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Pearson r (growth rate vs. CUE)")
    ax.set_ylim(-1.05, 1.05)
    ax.set_xticks(x)
    ax.set_title(
        f"Growth-CUE correlation across {n_max} substrates\n{subtitle}",
        fontsize=10,
        pad=8,
    )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            lw=0,
            markerfacecolor=LINE_C,
            markeredgecolor="white",
            markersize=7,
            label=f"p < {alpha:g}",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            lw=0,
            markerfacecolor="white",
            markeredgecolor=LINE_C,
            markersize=7,
            label=f"p $\\geq$ {alpha:g}",
        ),
        mpatch_proxy(),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="best")

    set_plot_style(ax)
    fig.tight_layout()
    fig.savefig(out_path / f"{filename}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_path / f"{filename}.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {filename}.png / .svg")


def mpatch_proxy():
    """Legend key for the bootstrap band."""
    from matplotlib.patches import Patch

    return Patch(
        facecolor=BAND_C,
        alpha=0.18,
        edgecolor="none",
        label="95% CI (bootstrap over substrates)",
    )


if __name__ == "__main__":
    main()
