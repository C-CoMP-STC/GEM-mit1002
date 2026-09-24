"""Tests for tools/release.py, the check that decides the release bump."""

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


class TestAgainstGit(unittest.TestCase):
    """End to end, against the last commit, so it needs git but no tags."""

    def test_committed_model_against_itself(self):
        result = release.check_release("HEAD")
        self.assertEqual(result.previous_path, MODEL_RELPATH)
        self.assertEqual(len(result.flips), 0)
        self.assertEqual(result.removed_reactions, [])


if __name__ == "__main__":
    unittest.main()
