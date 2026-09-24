"""Tests for tools/release.py, the check that decides the release bump."""

import os
import tempfile
import unittest

import pandas as pd

from tools import release
from tools.paths import MODEL_RELPATH


def _evaluation(calls, excluded=()):
    """A minimal stand-in for evaluate_phenotypes output: condition -> call."""
    return pd.DataFrame(
        {
            "condition": list(calls),
            "predicted": list(calls.values()),
            "growth": ["Yes"] * len(calls),
            "exclude_reason": [
                "not_representable" if c in excluded else "" for c in calls
            ],
        }
    )


class TestVersions(unittest.TestCase):
    def test_parse_accepts_old_v_prefix(self):
        self.assertEqual(release.parse_version("v3.1.0"), (3, 1, 0))
        self.assertEqual(release.parse_version("4.0.12"), (4, 0, 12))

    def test_parse_rejects_non_semver(self):
        for bad in ("4.0", "4.0.0-rc1", "four"):
            with self.assertRaises(ValueError):
                release.parse_version(bad)

    def test_bump(self):
        self.assertEqual(release.bump_version("3.1.0", "major"), "4.0.0")
        self.assertEqual(release.bump_version("3.1.4", "minor"), "3.2.0")
        self.assertEqual(release.bump_version("3.1.4", "patch"), "3.1.5")
        with self.assertRaises(ValueError):
            release.bump_version("3.1.0", "huge")

    def test_is_at_least(self):
        self.assertTrue(release.is_at_least("major", "minor"))
        self.assertTrue(release.is_at_least("minor", "minor"))
        self.assertFalse(release.is_at_least("patch", "minor"))

    def test_version_file_is_valid(self):
        release.read_version()  # raises if version.txt is malformed


class TestDiffGrowthCalls(unittest.TestCase):
    def test_no_change_no_flips(self):
        calls = {"a | x": "Yes", "b | y": "No"}
        flips = release.diff_growth_calls(_evaluation(calls), _evaluation(calls))
        self.assertEqual(len(flips), 0)

    def test_flip_in_either_direction(self):
        flips = release.diff_growth_calls(
            _evaluation({"a | x": "No", "b | y": "Yes"}),
            _evaluation({"a | x": "Yes", "b | y": "No"}),
        )
        self.assertEqual(sorted(flips["condition"]), ["a | x", "b | y"])
        row = flips.set_index("condition").loc["a | x"]
        self.assertEqual((row["previous"], row["current"]), ("No", "Yes"))
        self.assertTrue(row["agrees_now"])

    def test_invalid_solve_counts_as_a_flip(self):
        flips = release.diff_growth_calls(
            _evaluation({"a | x": "Yes"}), _evaluation({"a | x": None})
        )
        self.assertEqual(len(flips), 1)
        self.assertIsNone(flips.iloc[0]["current"])

    def test_excluded_rows_still_count(self):
        flips = release.diff_growth_calls(
            _evaluation({"a | x": "No"}, excluded={"a | x"}),
            _evaluation({"a | x": "Yes"}, excluded={"a | x"}),
        )
        self.assertEqual(len(flips), 1)
        self.assertTrue(flips.iloc[0]["excluded"])


def _check(**overrides):
    fields = {
        "against": "4.0.0",
        "previous_path": MODEL_RELPATH,
        "model_changed": False,
        "previous_id": "MIT1002_GEM",
        "current_id": "MIT1002_GEM",
        "flips": release.diff_growth_calls(_evaluation({}), _evaluation({})),
    }
    fields.update(overrides)
    return release.ReleaseCheck(**fields)


class TestRequiredBump(unittest.TestCase):
    def test_unchanged_model_is_patch(self):
        self.assertEqual(_check().required, "patch")

    def test_changed_model_without_flips_is_minor(self):
        self.assertEqual(_check(model_changed=True).required, "minor")

    def test_flip_is_major(self):
        flips = release.diff_growth_calls(
            _evaluation({"a | x": "No"}), _evaluation({"a | x": "Yes"})
        )
        self.assertEqual(_check(model_changed=True, flips=flips).required, "major")

    def test_moved_file_or_new_id_is_major(self):
        self.assertEqual(_check(previous_path="model.xml").required, "major")
        self.assertEqual(_check(previous_id="GEM_MIT1002").required, "major")

    def test_report_escapes_condition_pipes(self):
        flips = release.diff_growth_calls(
            _evaluation({"mbm | Glucose": "No"}), _evaluation({"mbm | Glucose": "Yes"})
        )
        report = release.format_report(_check(model_changed=True, flips=flips))
        self.assertIn("mbm \\| Glucose", report)


class TestChangelog(unittest.TestCase):
    def test_parse_conventional(self):
        self.assertEqual(
            release.parse_conventional("fix(model)!: drop rxn00001"),
            ("fix", "drop rxn00001", True),
        )
        self.assertEqual(
            release.parse_conventional("Untrack cache files"),
            ("other", "Untrack cache files", False),
        )

    def test_pr_title_from_merge_body(self):
        self.assertEqual(
            release.parse_pr_title(
                "fix: Ignore pseudoreactions\n", "fix/mass-imbalance"
            ),
            ("fix", "Ignore pseudoreactions", False),
        )

    def test_pr_title_made_from_branch_takes_kind_from_branch(self):
        self.assertEqual(
            release.parse_pr_title("Feat/gapfill pep", "feat/gapfill-pep"),
            ("feat", "gapfill pep", False),
        )

    def test_entry_splits_model_and_other_changes(self):
        changes = [
            release.Change("feat", "gapfill lysine", 425, touches_model=True),
            release.Change("docs", "fix typos", None, touches_model=False),
        ]
        flips = release.diff_growth_calls(
            _evaluation({"mbm | Lysine": "No"}), _evaluation({"mbm | Lysine": "Yes"})
        )
        check = _check(against="v3.1.0", model_changed=True, flips=flips)
        entry = release.changelog_entry("4.0.0", "2026-10-01", changes, check)
        self.assertTrue(entry.startswith("## 4.0.0 - 2026-10-01"))
        self.assertIn("Compared with 3.1.0: **major** release.", entry)
        model_part, other_part = entry.split("### Other changes")
        self.assertIn("gapfill lysine (#425)", model_part)
        self.assertIn("fix typos", other_part)
        self.assertIn("### Growth calls changed since 3.1.0", entry)

    def test_prepend_newest_first(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "CHANGELOG.md")
            release.prepend_changelog("## 1.0.0 - 2026-01-01\n", path)
            release.prepend_changelog("## 1.1.0 - 2026-02-01\n", path)
            with open(path, encoding="utf-8") as handle:
                text = handle.read()
        self.assertTrue(text.startswith("# Changelog"))
        self.assertLess(text.index("## 1.1.0"), text.index("## 1.0.0"))


class TestAgainstGit(unittest.TestCase):
    """End to end, against the last commit, so it needs git but no tags."""

    def test_committed_model_against_itself(self):
        result = release.check_release("HEAD")
        self.assertEqual(result.previous_path, MODEL_RELPATH)
        self.assertEqual(len(result.flips), 0)
        self.assertEqual(result.removed_reactions, [])

    def test_changes_since_parent_is_the_last_commit(self):
        changes = release.changes_since("HEAD~1")
        self.assertLessEqual(len(changes), 1)


if __name__ == "__main__":
    unittest.main()
