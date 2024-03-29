import psico.mcsalign
from pymol import cmd
import pytest


@pytest.mark.rdkit
def test_mcsalign__rdkit():
    cmd.reinitialize()
    cmd.fragment("his")
    cmd.fragment("trp")
    cmd.remove("hydro")
    psico.mcsalign.mcsalign("his", "trp", object="aln_r_0", method="rdkit")
    assert cmd.count_atoms("aln_r_0") == 20
    psico.mcsalign.mcsalign("his", "trp", object="aln_r_1", method="rdkit", exact=1)
    assert cmd.count_atoms("aln_r_1") == 18


def test_mcsalign__indigo():
    cmd.reinitialize()
    cmd.fragment("his")
    cmd.fragment("trp")
    cmd.remove("hydro")
    psico.mcsalign.mcsalign("his", "trp", object="aln_i_0", method="indigo")
    assert cmd.count_atoms("aln_i_0") == 18
    psico.mcsalign.mcsalign("his", "trp", object="aln_i_1", method="indigo", exact=1)
    assert cmd.count_atoms("aln_i_1") == 18
