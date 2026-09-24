"""Read a controlled vocabulary out of a markdown table in a README.

Several vocabularies in this repo are defined twice on purpose: once in code,
as the closed set a test enforces, and once in a README table, as the
explanation a contributor reads. The two have to agree, and the only way to
check that is to parse the table rather than search the document for the
words -- a substring search passes when a value happens to be mentioned
somewhere else, which is how ``duplicate`` and ``id_changed`` stayed
"documented" while nothing verified their table rows existed.

Used by ``test/test_deprecated.py`` and ``test/test_phenotype_data.py``.
"""

from __future__ import annotations

import os


def read_vocabulary_table(path: str, heading: str) -> dict[str, str]:
    """Return ``{value: description}`` from the first table under ``heading``.

    ``heading`` is the full markdown heading line, including its ``#`` marks,
    e.g. ``"### Controlled vocabulary for `reason`"``. Values are returned with
    surrounding backticks stripped, so a cell written ```` `duplicate` ````
    comes back as ``duplicate``.

    Raises ``AssertionError`` with a specific message when the file, the
    heading or the table is missing -- a silently empty result would make the
    tests that use this pass for the wrong reason.
    """
    assert os.path.exists(path), f"{path} is missing"

    lines = open(path, encoding="utf-8").read().splitlines()

    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == heading)
    except StopIteration:
        raise AssertionError(
            f"heading {heading!r} not found in {path}. If the heading was "
            f"renamed, update the constant in the test that reads it."
        ) from None

    # First table after the heading, stopping at the next heading so that a
    # section with no table fails loudly instead of picking up a later one.
    table_start = None
    for i in range(start + 1, len(lines)):
        stripped = lines[i].strip()
        if stripped.startswith("#"):
            break
        if stripped.startswith("|"):
            table_start = i
            break
    assert table_start is not None, (
        f"no markdown table under {heading!r} in {path}"
    )

    def cells(line: str) -> list[str]:
        return [c.strip() for c in line.strip().strip("|").split("|")]

    entries: dict[str, str] = {}
    # table_start is the header row; the row after it is the |---|---| rule.
    for line in lines[table_start + 2 :]:
        if not line.strip().startswith("|"):
            break
        row = cells(line)
        if not row or not row[0]:
            continue
        value = row[0].strip("`")
        # Re-join in case a description legitimately contains a pipe.
        description = "|".join(row[1:]).strip() if len(row) > 1 else ""
        entries[value] = description

    assert entries, f"table under {heading!r} in {path} has no rows"
    return entries
