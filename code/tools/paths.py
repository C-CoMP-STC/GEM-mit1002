"""Where things live in the MIT1002-GEM repository.

Every script, test and tool should get its locations from here instead of
counting parent directories from its own ``__file__``. Counting breaks as soon
as a file is moved, and it breaks quietly: the path still resolves, just to the
wrong folder. With the locations defined once, moving a folder means updating
this file and nothing else.

Usage::

    from tools.paths import DATA_DIR, MODEL_PATH

    model = cobra.io.read_sbml_model(MODEL_PATH)

A script outside ``tools/`` first has to make ``tools`` importable. It does
that by putting the folder that holds ``tools/`` on ``sys.path`` -- the one
place a script still counts parents -- and then imports from here::

    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from tools.paths import REPO_ROOT  # noqa: E402

All paths are absolute :class:`pathlib.Path` objects. Wrap them in ``str()``
for libraries that need a plain string, such as libSBML.
"""

from __future__ import annotations

from pathlib import Path

#: The ``code/`` folder, which holds ``tools/``, ``scripts/`` and ``test/``.
CODE_DIR = Path(__file__).resolve().parents[1]

#: Root of the git repository.
REPO_ROOT = CODE_DIR.parent

#: Name shared by the repository, the model files and (with ``_`` for ``-``)
#: the model ID, as standard-GEM requires.
MODEL_NAME = "MIT1002-GEM"

#: Folder holding the model in all its formats.
MODEL_DIR = REPO_ROOT / "model"

#: The SBML model, the source of truth for every other format.
MODEL_PATH = MODEL_DIR / f"{MODEL_NAME}.xml"

#: :data:`MODEL_PATH` relative to the repo root, as git and the GitHub API
#: see it.
MODEL_RELPATH = MODEL_PATH.relative_to(REPO_ROOT).as_posix()

#: Every repo-relative path the SBML model has had, newest first. It was
#: ``model.xml`` at the repo root until the standard-GEM rename, so anything
#: reading the model from an older commit or tag has to try both.
MODEL_RELPATH_HISTORY = (MODEL_RELPATH, "model.xml")

#: Version of the latest release, as standard-GEM requires.
VERSION_PATH = REPO_ROOT / "version.txt"

#: Experimental data, media provenance and derived tables.
DATA_DIR = REPO_ROOT / "data"

#: Genome sequences, gene calls and functional annotations of MIT1002.
GENOME_DIR = DATA_DIR / "genome"
