"""Schema tests for data/known_growth_phenotypes.tsv and its baseline.

Deliberately separate from ``test_growth.py``: none of this needs the model, so
it runs in milliseconds and still tells you which file is wrong when the slow
test fails. These guard the parts that other code assumes -- the exclusion
vocabulary, the uniqueness of the condition key, and the baseline pointing at
conditions that actually exist.

:class:`TestCountInterpretable` is the odd one out: a unit test of
``count_interpretable`` on a table built in the test, not a check on the
real file. It lives here because it needs neither the model nor the TSV.
"""

import csv
import os
import unittest

import pandas as pd

from tools.markdown_tables import read_vocabulary_table
from tools.phenotypes import (
    EXCLUSION_COLUMN,
    EXPECTED_MISMATCHES_TSV,
    EXCLUSION_REASONS,
    MISMATCH_COLUMNS,
    condition_key,
    count_interpretable,
    load_expected_mismatches,
    load_phenotypes,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_README = os.path.join(REPO_ROOT, "data", "README.md")

#: Heading of the table in data/README.md that defines the
#: ``exclude_reason`` vocabulary.
EXCLUSION_VOCAB_HEADING = "### `exclude_reason`"


class TestPhenotypeSchema(unittest.TestCase):
    def setUp(self):
        self.phenotypes = load_phenotypes()

    def test_growth_column_vocabulary(self):
        allowed = {"Yes", "No", "Unsure"}
        bad = sorted(set(self.phenotypes["growth"]) - allowed)
        self.assertFalse(
            bad,
            f"growth must be one of {sorted(allowed)}; found {bad}",
        )

    def test_exclude_reason_vocabulary(self):
        used = {r for r in self.phenotypes[EXCLUSION_COLUMN] if r}
        bad = sorted(used - set(EXCLUSION_REASONS))
        self.assertFalse(
            bad,
            f"{EXCLUSION_COLUMN} values not in the controlled vocabulary: {bad}. "
            f"Allowed: {sorted(EXCLUSION_REASONS)}. Adding a category means "
            f"updating tools/phenotypes.py and data/README.md together.",
        )

    def test_readme_table_matches_the_exclusion_vocabulary(self):
        """The vocabulary and the README have to be changed together.

        A set comparison against the table under ``### `exclude_reason```,
        not a substring search of the whole file: a value mentioned anywhere
        else must not count as documented, and a value dropped from
        ``EXCLUSION_REASONS`` must not be left in the table.
        """
        documented = set(read_vocabulary_table(DATA_README, EXCLUSION_VOCAB_HEADING))
        allowed = set(EXCLUSION_REASONS)
        self.assertEqual(
            allowed,
            documented,
            msg=(
                f"exclusion vocabulary and data/README.md disagree.\n"
                f"  allowed but undocumented: {sorted(allowed - documented)}\n"
                f"  documented but not allowed: {sorted(documented - allowed)}"
            ),
        )

    def test_every_documented_exclusion_reason_has_a_description(self):
        documented = read_vocabulary_table(DATA_README, EXCLUSION_VOCAB_HEADING)
        blank = sorted(v for v, desc in documented.items() if not desc)
        self.assertEqual(
            [], blank, msg=f"exclusion reasons listed with no explanation: {blank}"
        )

    def test_condition_keys_are_unique(self):
        """``minimal_media`` + ``c_source`` identifies a row.

        The expected-mismatch baseline is keyed on this pair, so a duplicate
        would make a baseline entry ambiguous.
        """
        counts = self.phenotypes["condition"].value_counts()
        duplicated = sorted(counts[counts > 1].index)
        self.assertFalse(
            duplicated,
            f"duplicate (minimal_media, c_source) keys: {duplicated}",
        )

    def test_met_ids_are_present_and_well_formed(self):
        bad = [
            (row["condition"], row["met_id"])
            for _, row in self.phenotypes.iterrows()
            if not row["met_id"]
            or not all(m.startswith("cpd") for m in row["met_id"])
        ]
        self.assertFalse(bad, f"rows with missing or malformed met_id: {bad}")


class TestCountInterpretable(unittest.TestCase):
    """Unit test for :func:`count_interpretable`, on a table built here.

    Deliberately not a check on ``known_growth_phenotypes.tsv``. There is no
    invariant to assert about the real file's count: it legitimately changes
    whenever an experiment comes in, and any value it takes is correct. The
    number quoted in the manuscript is pinned by the release the manuscript
    cites, not by this suite. What is worth testing is that the function
    implements the rule, which needs input whose answer is known by
    construction.
    """

    @staticmethod
    def _table(rows):
        """Build a phenotype-shaped frame from ``(growth, exclude_reason)`` pairs."""
        return pd.DataFrame(
            {
                "growth": [growth for growth, _ in rows],
                EXCLUSION_COLUMN: [reason for _, reason in rows],
            }
        )

    def test_counts_only_definite_unexcluded_rows(self):
        table = self._table(
            [
                ("Yes", ""),
                ("No", ""),
                ("Unsure", ""),
                ("Yes", "control_failed"),
                ("No", "id_uncertain"),
                ("Unsure", "conflicting_reports"),
            ]
        )
        self.assertEqual(count_interpretable(table), 2)

    def test_one_row_at_a_time(self):
        """Each condition for inclusion, in isolation, so a failure names itself."""
        cases = [
            (("Yes", ""), 1, "definite growth, not excluded"),
            (("No", ""), 1, "definite no-growth, not excluded"),
            (("Unsure", ""), 0, "no definite observation"),
            (("Yes", "control_failed"), 0, "excluded despite a definite observation"),
            (("Unsure", "control_failed"), 0, "held out for both reasons"),
            (("", ""), 0, "blank observation"),
        ]
        for row, expected, why in cases:
            with self.subTest(row=row, why=why):
                self.assertEqual(count_interpretable(self._table([row])), expected)

    def test_empty_table_is_zero(self):
        self.assertEqual(count_interpretable(self._table([])), 0)


class TestExpectedMismatches(unittest.TestCase):
    def setUp(self):
        self.phenotypes = load_phenotypes()
        self.expected = load_expected_mismatches()

    def test_baseline_has_the_expected_columns(self):
        """The baseline file must carry the columns its readers index by name.

        Checked against the file's own header row, not against the loaded
        frame. ``load_expected_mismatches`` backfills any absent column with
        empty strings so that an older baseline still loads, which means the
        frame always has every column and a test on the frame could never
        fail. The file is the thing a person edits or a script writes, so the
        file is the thing worth checking.

        A baseline missing ``minimal_media`` or ``c_source`` is the dangerous
        case: the loader fills them with "", every condition key becomes the
        same empty pair, and the baseline silently matches nothing -- so every
        accepted mismatch reads as new.
        """
        if not os.path.exists(EXPECTED_MISMATCHES_TSV):
            self.skipTest("no baseline file yet")
        with open(EXPECTED_MISMATCHES_TSV, newline="", encoding="utf-8") as handle:
            header = next(csv.reader(handle, delimiter="\t"), [])
        missing = [column for column in MISMATCH_COLUMNS if column not in header]
        self.assertFalse(
            missing,
            f"{os.path.relpath(EXPECTED_MISMATCHES_TSV, REPO_ROOT)} is missing "
            f"columns {missing}. load_expected_mismatches fills them in with "
            f"empty strings, so the file would load without complaint and "
            f"match nothing. Regenerate it with "
            f"scripts/update_phenotype_baseline.py.",
        )

    def test_baseline_rows_refer_to_real_conditions(self):
        """A baseline entry naming a condition that no longer exists is stale.

        This happens when a phenotype is renamed or removed and the baseline is
        not regenerated, and it would otherwise sit there forever excusing a
        mismatch that cannot occur.
        """
        known = set(self.phenotypes["condition"])
        orphaned = sorted(set(self.expected.get("condition", [])) - known)
        self.assertFalse(
            orphaned,
            f"baseline entries with no matching phenotype row: {orphaned}. "
            f"Regenerate with scripts/update_phenotype_baseline.py.",
        )

    def test_baseline_does_not_list_excluded_conditions(self):
        """Excluded rows are never scored, so they can never be a mismatch."""
        excluded = {
            condition_key(row["minimal_media"], row["c_source"])
            for _, row in self.phenotypes.iterrows()
            if row[EXCLUSION_COLUMN]
        }
        overlap = sorted(set(self.expected.get("condition", [])) & excluded)
        self.assertFalse(
            overlap,
            f"baseline lists conditions that are excluded from scoring, so the "
            f"entries are dead: {overlap}",
        )


if __name__ == "__main__":
    unittest.main()
