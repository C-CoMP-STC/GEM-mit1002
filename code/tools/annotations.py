"""Add MetaNetX cross-references to the annotations of the MIT1002 model.

Reaction IDs in the model are ModelSEED IDs with a compartment suffix
(``rxn00001_c0``). MetaNetX publishes a table, ``reac_xref.tsv``, that maps
IDs from other databases -- ModelSEED included -- onto MetaNetX reaction IDs
(``MNXR...``). This module reads that table and fills in the
``metanetx.reaction`` annotation for any reaction that does not have one.

Typical use, from the repo root::

    import cobra
    from tools.annotations import add_all_mnx_ids
    from tools.paths import MODEL_PATH

    model = cobra.io.read_sbml_model(MODEL_PATH)
    added = add_all_mnx_ids(model)
    print(f"added MetaNetX IDs to {len(added)} reactions")

Nothing is read when the module is imported. The cross-reference table is
loaded, and cached, the first time a lookup needs it.

``reac_xref.tsv`` is about 70 MB, so it is gitignored rather than committed.
Download the pinned release into the repo root with::

    curl -o reac_xref.tsv https://www.metanetx.org/ftp/4.4/reac_xref.tsv

The release is pinned (:data:`MNX_VERSION`) because MetaNetX IDs can be
merged or split between releases; a different version triggers a warning.
"""

from __future__ import annotations

import functools
import os
import types
import warnings
from collections import defaultdict
from collections.abc import Mapping

import cobra

from tools.paths import REPO_ROOT

# --------------------------------------------------------------------------
# MetaNetX release and file location
# --------------------------------------------------------------------------

#: MetaNetX release the annotations are taken from. Read from the ``#VERSION:``
#: line of the file and checked on load.
MNX_VERSION = "4.4"

#: Where to download the reaction cross-reference table for that release.
REAC_XREF_URL = f"https://www.metanetx.org/ftp/{MNX_VERSION}/reac_xref.tsv"

#: Default local copy of the table (gitignored).
REAC_XREF_PATH = os.path.join(REPO_ROOT, "reac_xref.tsv")

#: Annotation key (identifiers.org namespace) the MetaNetX IDs are stored under.
MNX_REACTION_NAMESPACE = "metanetx.reaction"

#: Placeholder MetaNetX uses for a source ID with no MetaNetX equivalent.
UNMAPPED_MNX_ID = "EMPTY"

#: A loaded cross-reference table: source ID (no prefix) -> MetaNetX IDs.
XrefTable = Mapping[str, tuple[str, ...]]


# --------------------------------------------------------------------------
# Reading the cross-reference table
# --------------------------------------------------------------------------


@functools.cache
def load_reac_xref(
    source: str = "seed.reaction", path: str = REAC_XREF_PATH
) -> XrefTable:
    """Read the MetaNetX reaction cross-references for one source database.

    Only rows for ``source`` are kept, which keeps memory small even though
    the file covers every database MetaNetX reconciles. The result is cached
    per ``(source, path)``, so repeated lookups do not re-read the file.

    Args:
        source: Database prefix as written in the first column of
            ``reac_xref.tsv``, e.g. ``"seed.reaction"`` or ``"bigg.reaction"``.
        path: Location of ``reac_xref.tsv``.

    Returns:
        A read-only mapping from source IDs, with the ``source:`` prefix
        removed (``"rxn00001"``), to the MetaNetX reaction IDs they map to.
        Source IDs that MetaNetX marks as unmapped (``EMPTY``) are omitted.

    Raises:
        FileNotFoundError: If ``path`` does not exist. The message says how to
            download the file.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"MetaNetX cross-reference table not found at {path}. "
            f"Download it with:\n    curl -o {path} {REAC_XREF_URL}"
        )

    prefix = source + ":"
    table: defaultdict[str, list[str]] = defaultdict(list)
    version: str | None = None

    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                if line.startswith("#VERSION:"):
                    version = line.split(":", 1)[1].strip()
                continue
            if not line.startswith(prefix):
                continue
            source_ref, mnx_id, _description = line.rstrip("\n").split("\t", 2)
            if mnx_id == UNMAPPED_MNX_ID:
                continue
            mnx_ids = table[source_ref[len(prefix) :]]
            if mnx_id not in mnx_ids:
                mnx_ids.append(mnx_id)

    if version != MNX_VERSION:
        warnings.warn(
            f"{path} is MetaNetX version {version!r}, but annotations are "
            f"pinned to {MNX_VERSION!r}. IDs may differ between releases.",
            stacklevel=2,
        )

    # Read-only, because this object is shared through the cache.
    return types.MappingProxyType({k: tuple(v) for k, v in table.items()})


# --------------------------------------------------------------------------
# Lookups
# --------------------------------------------------------------------------


def strip_compartment(identifier: str) -> str:
    """Remove the compartment suffix from a model ID.

    Args:
        identifier: A model ID such as ``"rxn00001_c0"``.

    Returns:
        The ID without its final ``_<compartment>`` part (``"rxn00001"``).
        An ID with no underscore is returned unchanged.
    """
    base, separator, _compartment = identifier.rpartition("_")
    return base if separator else identifier


def get_mnx_id(
    identifier: str,
    source: str = "seed.reaction",
    xref: XrefTable | None = None,
) -> tuple[str, ...]:
    """Look up the MetaNetX reaction IDs for one source-database ID.

    Args:
        identifier: ID in the source database, without compartment suffix or
            ``source:`` prefix, e.g. ``"rxn00001"``.
        source: Source database prefix, as for :func:`load_reac_xref`.
        xref: A table from :func:`load_reac_xref`. Loaded (and cached) from
            the default path if not given.

    Returns:
        The matching MetaNetX IDs. Empty if there is no match; more than one
        if MetaNetX maps the source ID to several reactions.
    """
    if xref is None:
        xref = load_reac_xref(source)
    return xref.get(identifier, ())


# --------------------------------------------------------------------------
# Annotating a model
# --------------------------------------------------------------------------


def add_all_mnx_ids(
    model: cobra.Model,
    source: str = "seed.reaction",
    xref: XrefTable | None = None,
) -> dict[str, tuple[str, ...]]:
    """Add MetaNetX IDs to reactions that do not already have one.

    Existing ``metanetx.reaction`` annotations are never changed, and boundary
    (exchange, demand, sink) reactions are skipped because they have no
    database equivalent. The model is modified in place; writing it back to
    disk is left to the caller.

    A single match is stored as a string. Several matches are all stored, as
    a list, so that the ambiguity stays visible in the model rather than one
    ID being picked arbitrarily.

    Args:
        model: The model to annotate.
        source: Database the reaction IDs come from, as for
            :func:`load_reac_xref`.
        xref: A table from :func:`load_reac_xref`. Loaded (and cached) from
            the default path if not given.

    Returns:
        The reactions that were annotated, as reaction ID -> MetaNetX IDs
        added. Useful for reporting the change in a commit or PR.
    """
    if xref is None:
        xref = load_reac_xref(source)

    added: dict[str, tuple[str, ...]] = {}
    for reaction in model.reactions:
        if reaction.boundary or MNX_REACTION_NAMESPACE in reaction.annotation:
            continue
        mnx_ids = get_mnx_id(strip_compartment(reaction.id), source, xref)
        if not mnx_ids:
            continue
        reaction.annotation[MNX_REACTION_NAMESPACE] = (
            mnx_ids[0] if len(mnx_ids) == 1 else list(mnx_ids)
        )
        added[reaction.id] = mnx_ids
    return added
