"""Tests for tools/annotations.py.

These use a small hand-written cross-reference file, so they do not need the
70 MB ``reac_xref.tsv`` to be downloaded.
"""

import warnings

import cobra
import pytest

from tools.annotations import (
    MNX_REACTION_NAMESPACE,
    add_all_mnx_ids,
    get_mnx_id,
    load_reac_xref,
    strip_compartment,
)

XREF_TEXT = """\
### MetaNetX/MNXref reconciliation ###
#VERSION:   4.4
seed.reaction:rxn00001\tMNXR1\tone match
seed.reaction:rxn00002\tMNXR2a\ttwo matches
seed.reaction:rxn00002\tMNXR2b\ttwo matches
seed.reaction:rxn00003\tEMPTY\tunmapped
seedR:rxn00004\tMNXR4\tother prefix only
bigg.reaction:PPA\tMNXR1\tanother database
"""


@pytest.fixture
def xref_path(tmp_path):
    path = tmp_path / "reac_xref.tsv"
    path.write_text(XREF_TEXT)
    return str(path)


def test_load_keeps_only_requested_source(xref_path):
    table = load_reac_xref("seed.reaction", xref_path)
    assert dict(table) == {"rxn00001": ("MNXR1",), "rxn00002": ("MNXR2a", "MNXR2b")}


def test_load_warns_on_other_version(tmp_path):
    path = tmp_path / "reac_xref.tsv"
    path.write_text(XREF_TEXT.replace("4.4", "4.3"))
    with pytest.warns(UserWarning, match="4.3"):
        load_reac_xref("seed.reaction", str(path))


def test_load_missing_file_says_how_to_download(tmp_path):
    with pytest.raises(FileNotFoundError, match="curl"):
        load_reac_xref("seed.reaction", str(tmp_path / "absent.tsv"))


def test_strip_compartment():
    assert strip_compartment("rxn00001_c0") == "rxn00001"
    assert strip_compartment("rxn00001") == "rxn00001"


def test_get_mnx_id_no_match_is_empty(xref_path):
    table = load_reac_xref("seed.reaction", xref_path)
    assert get_mnx_id("rxn99999", xref=table) == ()


def test_add_all_mnx_ids(xref_path):
    table = load_reac_xref("seed.reaction", xref_path)
    model = cobra.Model("test")
    single, multi, unmapped, annotated = (
        cobra.Reaction(i)
        for i in ("rxn00001_c0", "rxn00002_c0", "rxn00003_c0", "rxn00005_c0")
    )
    annotated.annotation[MNX_REACTION_NAMESPACE] = "MNXR_existing"
    met = cobra.Metabolite("cpd00001_c0", compartment="c0")
    for rxn in (single, multi, unmapped, annotated):
        rxn.add_metabolites(
            {met: -1, cobra.Metabolite(rxn.id + "_p", compartment="c0"): 1}
        )
    model.add_reactions([single, multi, unmapped, annotated])

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        added = add_all_mnx_ids(model, xref=table)

    assert added == {"rxn00001_c0": ("MNXR1",), "rxn00002_c0": ("MNXR2a", "MNXR2b")}
    assert single.annotation[MNX_REACTION_NAMESPACE] == "MNXR1"
    assert multi.annotation[MNX_REACTION_NAMESPACE] == ["MNXR2a", "MNXR2b"]
    assert MNX_REACTION_NAMESPACE not in unmapped.annotation
    assert annotated.annotation[MNX_REACTION_NAMESPACE] == "MNXR_existing"
