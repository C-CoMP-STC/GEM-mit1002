"""Check that tools/paths.py points at folders that exist.

If a folder is moved and ``tools/paths.py`` is not updated, every script
fails with a confusing "file not found". These tests fail first, and say why.
"""

import unittest

from tools import paths


class TestPaths(unittest.TestCase):
    def test_code_dir_holds_the_code(self):
        for folder in ("tools", "scripts", "test"):
            self.assertTrue(
                (paths.CODE_DIR / folder).is_dir(),
                f"{folder}/ is not in CODE_DIR ({paths.CODE_DIR})",
            )

    def test_repo_root_is_the_repository_root(self):
        # .github/ rather than .git/, so this also passes in a downloaded
        # release, which has no .git folder.
        self.assertTrue(
            (paths.REPO_ROOT / ".github").is_dir(),
            f"REPO_ROOT ({paths.REPO_ROOT}) is not the root of the repository",
        )

    def test_model_and_data_exist(self):
        self.assertTrue(paths.MODEL_PATH.is_file(), f"missing {paths.MODEL_PATH}")
        self.assertTrue(paths.DATA_DIR.is_dir(), f"missing {paths.DATA_DIR}")
        self.assertTrue(paths.GENOME_DIR.is_dir(), f"missing {paths.GENOME_DIR}")

    def test_model_relpath(self):
        self.assertEqual(paths.MODEL_RELPATH, "model/MIT1002-GEM.xml")


if __name__ == "__main__":
    unittest.main()
