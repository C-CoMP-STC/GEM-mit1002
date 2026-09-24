"""Consistency checks on the committed curation time series.

``code/curation_process/phenotype_confusion_over_time.csv`` is the file figure 2 of
the manuscript is read off: panel B plots its match counts and panel A quotes
two of its rows as confusion matrices. Nothing recomputes it in CI -- producing
it needs the GitHub API and a solve per PR -- so what is checked here is that
the committed file is internally coherent and was produced by one scoring
definition.

That is the failure this catches. The numbers in the figure originally could not
be reconciled with the phenotype table by hand, and the reason was a file
holding rows scored under two different definitions of "match", with no record
of which was which.

No model is loaded, so this runs in milliseconds.
"""

import os
import unittest

import pandas as pd

from tools.paths import CODE_DIR
from tools.phenotypes import count_interpretable, load_phenotypes

CONFUSION_CSV = os.path.join(
    CODE_DIR, "curation_process", "phenotype_confusion_over_time.csv"
)
SUMMARY_CSV = os.path.join(CODE_DIR, "curation_process", "growth_match_summary.csv")

#: Kept in sync with ``run_tests_on_prs.SCORING_VERSION`` by
#: :meth:`TestConfusionTimeline.test_scoring_version_matches_the_script`, rather
#: than imported, so that this file does not need cobra installed.
EXPECTED_SCORING_VERSION = 3

#: PRs whose model file cannot be read by cobra, so they are recorded as
#: ERROR on every run. These are real gaps in the history, not failed
#: downloads, and re-running the script will not fill them. Each entry says why,
#: so that the gap in figure 2B can be explained.
#:
#: 285-289: commit 3a146b7 (in #285) added the 1,2-ethanediol transporter and
#: EX_cpd00992_e0 but not the species M_cpd00992_e0, so reading the SBML fails
#: with KeyError 'cpd00992_e0'. #286, #288 and #289 branched from dev after #285
#: and inherited it. Fixed by e65558b in #290.
KNOWN_INVALID_MODELS = {
    285: "M_cpd00992_e0 referenced but not defined; fixed in #290",
    286: "M_cpd00992_e0 referenced but not defined (inherited from #285); fixed in #290",
    288: "M_cpd00992_e0 referenced but not defined (inherited from #285); fixed in #290",
    289: "M_cpd00992_e0 referenced but not defined (inherited from #285); fixed in #290",
}


def _load(path):
    table = pd.read_csv(path)
    # Rows whose evaluation failed are recorded as ERROR on purpose; they carry
    # no counts to check, and dropping them here is not hiding anything because
    # test_no_errored_rows reports them separately.
    table = table[table["Matches"].astype(str) != "ERROR"].copy()
    # A single ERROR row makes pandas read every column after it as strings,
    # and filtering the row out does not change the dtype back. Without this,
    # the arithmetic below concatenates ("19" + "33") instead of adding.
    numeric = table.columns.difference(["Date Opened", "Date Merged"])
    table[numeric] = table[numeric].apply(pd.to_numeric)
    return table


#: What to tell someone when a timeline file is missing. These files are
#: committed and figure 2 is read off them, so a missing file means a broken
#: path or a bad commit, never "not generated yet" -- and it must fail, not skip.
MISSING_FILE_HELP = (
    "{name} is missing from code/curation_process/. It is committed and "
    "figure 2 depends on it, so either it was deleted or tools.paths no longer "
    "points at it. Restore it from git, or regenerate it by running "
    "code/curation_process/run_tests_on_prs.py (needs `gh` authenticated)."
)


class TestTimelineFilesExist(unittest.TestCase):
    """Fail loudly, with one clear message, if a timeline file is missing."""

    def test_confusion_csv_exists(self):
        self.assertTrue(
            os.path.isfile(CONFUSION_CSV),
            MISSING_FILE_HELP.format(name=os.path.basename(CONFUSION_CSV)),
        )

    def test_summary_csv_exists(self):
        self.assertTrue(
            os.path.isfile(SUMMARY_CSV),
            MISSING_FILE_HELP.format(name=os.path.basename(SUMMARY_CSV)),
        )


class TestConfusionTimeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.path.isfile(CONFUSION_CSV):
            raise FileNotFoundError(
                MISSING_FILE_HELP.format(name=os.path.basename(CONFUSION_CSV))
            )
        cls.raw = pd.read_csv(CONFUSION_CSV)
        cls.data = _load(CONFUSION_CSV)
        cls.phenotypes = load_phenotypes()

    def test_one_scoring_version_throughout(self):
        versions = sorted(self.raw["Scoring Version"].dropna().unique())
        self.assertEqual(
            versions,
            [EXPECTED_SCORING_VERSION],
            "the time series mixes scoring definitions (or predates the "
            "version stamp), so its points are not comparable with each other. "
            "Re-run code/curation_process/run_tests_on_prs.py.",
        )

    def test_scoring_version_matches_the_script(self):
        """Catch the stamp and the constant drifting apart."""
        script = os.path.join(
            CODE_DIR, "curation_process", "run_tests_on_prs.py"
        )
        with open(script, encoding="utf-8") as handle:
            source = handle.read()
        self.assertIn(
            f"SCORING_VERSION = {EXPECTED_SCORING_VERSION}",
            source,
            "run_tests_on_prs.py declares a different SCORING_VERSION than this "
            "test expects; bump EXPECTED_SCORING_VERSION and regenerate the CSV",
        )

    def test_confusion_matrix_sums_to_scored(self):
        cells = (
            self.data["True Positive"]
            + self.data["True Negative"]
            + self.data["False Positive"]
            + self.data["False Negative"]
        )
        bad = self.data.loc[cells != self.data["Scored"], "PR Number"].tolist()
        self.assertFalse(bad, f"TP+TN+FP+FN != Scored for PR(s) {bad}")

    def test_matches_equal_the_concordant_cells(self):
        concordant = self.data["True Positive"] + self.data["True Negative"]
        bad = self.data.loc[
            concordant != self.data["Matches"], "PR Number"
        ].tolist()
        self.assertFalse(bad, f"Matches != TP+TN for PR(s) {bad}")

    def test_all_categories_sum_to_the_condition_count(self):
        """Every condition lands in exactly one bucket, for every PR."""
        total = (
            self.data["Scored"]
            + self.data["Unsure"]
            + self.data["Excluded"]
            + self.data["Invalid Solve"]
        )
        bad = self.data.loc[total != self.data["Conditions"], "PR Number"].tolist()
        self.assertFalse(
            bad,
            f"the category counts do not add up to Conditions for PR(s) {bad}; "
            f"a condition is being double counted or dropped",
        )

    def test_no_uptake_route_is_within_the_negative_predictions(self):
        """It is a subset of Scored, not another bucket.

        Quoting it next to the confusion matrix only works if it fits inside
        TN + FN. A count larger than that means unsure or excluded rows leaked
        into it.
        """
        negatives = self.data["True Negative"] + self.data["False Negative"]
        bad = self.data.loc[
            self.data["No Uptake Route"] > negatives, "PR Number"
        ].tolist()
        self.assertFalse(
            bad, f"No Uptake Route exceeds TN+FN for PR(s) {bad}"
        )

    def test_condition_count_matches_the_phenotype_table(self):
        """A stale series scored against a different version of the TSV.

        Rows are added to ``known_growth_phenotypes.tsv`` as experiments come
        in. Because the series is cached per PR, an old row can have been
        scored against a smaller table -- which makes the line's early points
        incomparable with its late ones for a reason that has nothing to do
        with the model.
        """
        expected = len(self.phenotypes)
        bad = self.data.loc[
            self.data["Conditions"] != expected, "PR Number"
        ].tolist()
        self.assertFalse(
            bad,
            f"PR(s) {bad} were scored against a different number of conditions "
            f"than the {expected} now in data/known_growth_phenotypes.tsv. "
            f"Re-run run_tests_on_prs.py with FORCE_RERUN = True.",
        )

    def test_denominator_is_constant_and_correct(self):
        expected = count_interpretable(self.phenotypes)
        found = sorted(self.data["Interpretable"].unique())
        self.assertEqual(
            found,
            [expected],
            f"the match denominator should be {expected} for every PR (the "
            f"conditions with a definite Yes/No and no exclusion reason); "
            f"found {found}",
        )

    def test_matches_within_range(self):
        bad = self.data.loc[
            (self.data["Matches"] < 0)
            | (self.data["Matches"] > self.data["Interpretable"]),
            "PR Number",
        ].tolist()
        self.assertFalse(bad, f"Matches outside 0..Interpretable for PR(s) {bad}")

    def test_pr_numbers_are_unique(self):
        counts = self.raw["PR Number"].value_counts()
        self.assertFalse(
            sorted(counts[counts > 1].index),
            f"duplicate PR rows: {sorted(counts[counts > 1].index)}",
        )

    def test_no_errored_rows(self):
        """An ERROR row is a gap in the figure, not a result.

        Not fatal to the analysis, but it should be visible rather than sitting
        in the file unnoticed -- re-running the script retries them.

        PRs in :data:`KNOWN_INVALID_MODELS` are allowed: their models cannot be
        read, so ERROR is the correct record for them.
        """
        errored = self.raw.loc[
            self.raw["Matches"].astype(str) == "ERROR", "PR Number"
        ].tolist()
        unexplained = [pr for pr in errored if pr not in KNOWN_INVALID_MODELS]
        self.assertFalse(
            unexplained,
            f"PR(s) {unexplained} failed to evaluate and are missing from the "
            f"series; re-run code/curation_process/run_tests_on_prs.py. If the "
            f"model at that PR cannot be read, add it to KNOWN_INVALID_MODELS "
            f"with the reason.",
        )

    def test_known_invalid_models_are_still_errored(self):
        """Keep KNOWN_INVALID_MODELS from going stale.

        If one of these PRs now has scores, either the model became readable
        (e.g. a cobra update) or the history was rewritten. Either way the
        entry no longer describes the file and should be removed.
        """
        rows = self.raw.set_index("PR Number")
        missing = [pr for pr in KNOWN_INVALID_MODELS if pr not in rows.index]
        self.assertFalse(
            missing,
            f"KNOWN_INVALID_MODELS lists PR(s) {missing} that are not in the "
            f"time series at all",
        )
        scored = [
            pr
            for pr in KNOWN_INVALID_MODELS
            if str(rows.loc[pr, "Matches"]) != "ERROR"
        ]
        self.assertFalse(
            scored,
            f"PR(s) {scored} are listed in KNOWN_INVALID_MODELS but now have "
            f"scores; remove them from the list",
        )

    def test_summary_view_agrees_with_the_full_record(self):
        """``growth_match_summary.csv`` is a view, so it must not disagree."""
        if not os.path.isfile(SUMMARY_CSV):
            self.fail(MISSING_FILE_HELP.format(name=os.path.basename(SUMMARY_CSV)))
        summary = _load(SUMMARY_CSV).set_index("PR Number")
        full = self.data.set_index("PR Number")
        shared = summary.index.intersection(full.index)
        self.assertEqual(
            len(shared),
            len(full),
            "the two files cover different sets of PRs; regenerate both by "
            "running run_tests_on_prs.py",
        )
        for column, source in (("Matches", "Matches"), ("Total", "Interpretable")):
            mismatched = shared[
                summary.loc[shared, column].values != full.loc[shared, source].values
            ].tolist()
            self.assertFalse(
                mismatched,
                f"{column} in growth_match_summary.csv disagrees with "
                f"{source} in phenotype_confusion_over_time.csv for PR(s) "
                f"{mismatched}",
            )


if __name__ == "__main__":
    unittest.main()
