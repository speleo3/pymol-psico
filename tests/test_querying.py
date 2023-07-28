import psico.querying
from pymol import cmd, CmdException
import pytest
from pytest import approx


def test_get_color():
    cmd.reinitialize()
    cmd.pseudoatom(color="0xFF0000")
    cmd.pseudoatom(color="0xFFFF00")  # most frequent
    cmd.pseudoatom(color="0x00FF00")  # middle
    cmd.pseudoatom(color="0xFFFF00")  # most frequent
    cmd.pseudoatom(color="0x0000FF")
    assert psico.querying.get_color("*", which=0) == "0xff0000"
    assert psico.querying.get_color("*", which=1) == "0x00ff00"
    assert psico.querying.get_color("*", which=2) == "0xffff00"
    assert psico.querying.get_color("*", mode=0) == "0xff0000"
    assert psico.querying.get_color("*", mode=1) == (1.0, 0.0, 0.0)
    assert psico.querying.get_color("*", mode=2) == "#ff0000"
    assert psico.querying.get_color("none") == "gray"


def test_iterate_to_list():
    cmd.reinitialize()
    cmd.fragment("gly")
    r = psico.querying.iterate_to_list("not hydro", "name")
    assert r == ["N", "CA", "C", "O"]


def test_iterate_to_list__space():
    cmd.reinitialize()
    cmd.pseudoatom()
    space = {"myspacevar": 123}
    r = psico.querying.iterate_to_list("all", "myspacevar", space=space)
    assert r == [123]


def test_iterate_state_to_list():
    cmd.reinitialize()
    cmd.fragment("gly")
    r = psico.querying.iterate_state_to_list(-1, "not hydro", "(x,y,z)")
    assert len(r) == 4
    assert r[0] == approx([-1.1945616, 0.201105937, -0.2060566])
    assert r[1] == approx([0.23043768, 0.318106502, -0.5020566])
    assert r[2] == approx([1.05943620, -0.38989559, 0.54194343])
    assert r[3] == approx([0.54543632, -0.97489464, 1.49894345])


def test_iterate_to_sele():
    cmd.reinitialize()
    cmd.fab("AGAG", "m1")
    sele = psico.querying.iterate_to_sele("m1", "resn == 'GLY'")
    assert cmd.count_atoms(sele) == cmd.count_atoms("resn GLY")


def test_csp__pymol2():
    from pymol2 import PyMOL
    with PyMOL() as p1, PyMOL() as p2:
        p1.cmd.fab("A/1/ ARKA B/1/ GERD", "m1")
        p2.cmd.fab("EE", "m1")
        p2.cmd.fab("DD", "m2")
        assert psico.querying.csp("m1", _self=p1.cmd) == -2
        assert psico.querying.csp("m1", "m2", _self=p2.cmd) == 4


def test_csp():
    cmd.reinitialize()
    cmd.fab("A/1/ ARKA B/1/ GERD", "m1")
    cmd.fab("EE", "m2")
    cmd.fab("DD", "m3")
    assert psico.querying.csp("m1") == -2
    assert psico.querying.csp("m2", "m3") == 4


def test_extinction_coefficient():
    cmd.reinitialize()
    eps_reduced = 2 * 5500 + 3 * 1490
    eps_ss_1 = eps_reduced + 1 * 125
    eps_ss_2 = eps_reduced + 2 * 125
    cmd.fab("AWWYYYCCCC", "m1")
    eps, A_280 = psico.querying.extinction_coefficient("m1")
    assert eps == eps_reduced
    assert A_280 == approx(eps_reduced / 1345.58908)
    cmd.remove("hydro")
    eps, A_280 = psico.querying.extinction_coefficient("m1")
    assert eps == eps_reduced
    assert A_280 == approx(eps_reduced / 1347.60496)
    cmd.bond("7/SG", "8/SG")
    eps, A_280 = psico.querying.extinction_coefficient("m1")
    assert eps == eps_ss_1
    assert A_280 == approx(eps_ss_1 / 1345.58908)
    cmd.bond("9/SG", "10/SG")
    eps, A_280 = psico.querying.extinction_coefficient("m1")
    assert eps == eps_ss_2
    assert A_280 == approx(eps_ss_2 / 1343.5732)


def test_isoelectric_point():
    cmd.reinitialize()
    cmd.fab("AWWYYYCCCC", "m1")
    cmd.fab("DDDKKEEEHHHIIKL", "m2")
    pi = psico.querying.isoelectric_point()
    assert pi == approx(5.31, abs=1e-2)


def test_get_segis():
    cmd.reinitialize()
    cmd.pseudoatom("m1", segi="seg1")
    cmd.pseudoatom("m2", segi="seg2")
    assert psico.querying.get_segis() == {"seg1", "seg2"}


def test_centerofmass():
    cmd.reinitialize()
    cmd.fragment("gly")
    cmd.alter("all", "q=index/7.0")
    r = psico.querying.centerofmass()
    assert r == approx([0.384, -0.389, 0.578], abs=1e-3)


def test_gyradius():
    cmd.reinitialize()
    cmd.fragment("gly")
    cmd.alter("all", "q=index/7.0")
    r = psico.querying.gyradius()
    assert r == approx(1.28512885)


def test_get_alignment_coords():
    cmd.reinitialize()
    cmd.fragment("gly")
    cmd.fragment("val")
    cmd.fragment("ala")
    cmd.align("gly", "ala & !hydro", object="aln")
    cmd.align("val", "ala & !hydro", object="aln")
    assert cmd.count_atoms("aln & gly") == 4
    assert cmd.count_atoms("aln & val") == 5
    r = psico.querying.get_alignment_coords("aln")
    assert set(r) == {"ala", "val", "gly"}
    assert all(len(c) == 4 for c in r.values())
    assert r["ala"][0] == approx([-0.6769, -1.2303, -0.4905], abs=1e-4)
    assert r["ala"][1] == approx([-0.0009, 0.0637, -0.4905], abs=1e-4)
    assert r["gly"][0] == approx([-0.67673, -1.23049, -0.49175], abs=1e-4)
    assert r["gly"][3] == approx([2.02836, -1.22582, -0.49986], abs=1e-4)


def test_get_sasa():
    cmd.reinitialize()
    cmd.fragment("gly")
    r = psico.querying.get_sasa("gly", dot_density=4)
    assert r == approx(200.888, abs=1e-3)
    assert cmd.get_setting_int("dot_solvent") == 0
    assert cmd.get_setting_int("dot_density") == 2
    assert cmd.get_setting_int("dot_density", "gly") == 2


def test_get_object_state():
    cmd.reinitialize()
    cmd.fragment("gly")
    cmd.fragment("ala")
    with pytest.raises(CmdException, match="No object"):
        psico.querying.get_object_name("none")
    with pytest.raises(CmdException, match="more than one object"):
        psico.querying.get_object_name("all", strict=1)
    assert psico.querying.get_object_name("all") == "gly"
    assert psico.querying.get_object_name("name CB") == "ala"


def test_get_selection_state():
    cmd.reinitialize()
    cmd.fragment("cys")
    cmd.fragment("ser")
    cmd.create("cys", "cys", 1, 2)
    cmd.create("cys", "cys", 1, 3)
    assert psico.querying.get_selection_state("name SG") == 1
    assert psico.querying.get_selection_state("name OG") == 1
    cmd.frame(3)
    assert psico.querying.get_selection_state("name SG") == 3
    assert psico.querying.get_selection_state("name OG") == 1
    cmd.set("state", 2, "cys")
    cmd.set("static_singletons", 0)
    assert psico.querying.get_selection_state("name SG") == 2
    with pytest.raises(CmdException, match="Invalid state"):
        psico.querying.get_selection_state("name OG")


def test_get_raw_distances():
    cmd.reinitialize()
    cmd.fragment("gly")
    cmd.distance("d1", "gly & elem C", "gly & elem N")
    d = sorted(psico.querying.get_raw_distances("d1"))
    assert len(d) == 2
    assert d[0] == (('gly', 2), ('gly', 1), approx(1.4601, abs=1e-3))
    assert d[1] == (('gly', 3), ('gly', 1), approx(2.4473, abs=1e-3))


def test_shortest_distance():
    cmd.reinitialize()
    cmd.fragment(name="ala", origin=0)
    cmd.fragment(name="his", origin=0)
    dist, ala_atom, his_atom = psico.querying.shortest_distance("ala", "his")
    assert dist == approx(7.4, rel=1e-2)
    assert ala_atom == "/ala///ALA`2/2HB"
    assert his_atom == "/his///HIS`2/O"
