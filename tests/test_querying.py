import psico.querying
from pymol import cmd
from pytest import approx


def test_iterate_to_list():
    cmd.reinitialize()
    cmd.fragment("gly")
    r = psico.querying.iterate_to_list("not hydro", "name")
    assert r == ["N", "CA", "C", "O"]


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


#def centerofmass(selection='(all)', state=-1, quiet=1):
#def gyradius(selection='(all)', state=-1, quiet=1):
#def get_alignment_coords(name, active_only=0, state=-1, quiet=0):
#def get_sasa(selection, state=-1, dot_density=5, quiet=1):
#def get_sasa_ball(selection, state=-1, quiet=1):
#def get_sasa_mmtk(selection, state=-1, hydrogens='auto', quiet=1):
#def get_raw_distances(names='', state=1, selection='all', quiet=1):
#def get_color(selection, which=0, mode=0):
#def get_object_name(selection, strict=0):
#def get_object_state(name):
#def get_selection_state(selection):
#def get_ensemble_coords(selection):

def test_shortest_distance():
    cmd.reinitialize()
    cmd.fragment(name="ala", origin=0)
    cmd.fragment(name="his", origin=0)
    dist, ala_atom, his_atom = psico.querying.shortest_distance("ala", "his")
    assert dist == approx(7.4, rel=1e-2)
    assert ala_atom == "/ala///ALA`2/2HB"
    assert his_atom == "/his///HIS`2/O"
