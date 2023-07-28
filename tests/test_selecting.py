import psico.selecting
from pymol import cmd


def test_select_pepseq():
    cmd.reinitialize()
    cmd.fab("ACDEFGHCEEA")
    psico.selecting.select_pepseq("CDE")
    assert cmd.count_atoms("sele") == 38
    psico.selecting.select_pepseq("C[DE]E")
    assert cmd.count_atoms("sele") == 79
    psico.selecting.select_pepseq("CE*", name="s1")
    assert cmd.count_atoms("s1") == 52
    psico.selecting.select_pepseq("CE+", name="s1")
    assert cmd.count_atoms("s1") == 41
    psico.selecting.select_pepseq("X.Z", one_letter={"PHE": "X", "HIS": "Z"})
    assert cmd.count_atoms("sele") == 44


def test_diff():
    cmd.reinitialize()
    cmd.fab("AGGG", "m1")
    cmd.fab("GGGS", "m2")
    r = psico.selecting.diff("m1", "m2")
    assert r == "diff01"
    assert cmd.count_atoms(r) == 17
    r = psico.selecting.diff("m1", "m2", operator="like")
    assert r == "diff02"
    assert cmd.count_atoms(r) == 5
    r = psico.selecting.symdiff("m1", "m2", name="s2")
    assert r == "s2"
    assert cmd.count_atoms(r) == 35


def test_collapse_resi():
    cmd.reinitialize()
    cmd.fab("ACDEFGHI", "m1")
    r = psico.selecting.collapse_resi("resn ALA+CYS+ASP+GLY")
    assert r == "/m1///1-3+6"
    r = psico.selecting.collapse_resi("name CG")
    assert r == "/m1///3-5+7"


def test_select_distances():
    cmd.reinitialize()
    cmd.fragment("gly")
    cmd.distance("d1", "gly & elem C", "gly & elem N")
    psico.selecting.select_distances("d1", "s1")
    assert cmd.count_atoms("s1") == 3


def test_select_range():
    cmd.reinitialize()
    cmd.fab("GASGAGS", "m1")
    cmd.select("s1", "resn ALA")
    psico.selecting.select_range()
    assert 4 == cmd.count_atoms("s1 & guide")
    psico.selecting.select_range("s2", "resn SER")
    assert 5 == cmd.count_atoms("s2 & guide")
