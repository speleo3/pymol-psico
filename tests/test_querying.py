import psico.querying
from pymol import cmd


def test_iterate_to_list():
    cmd.reinitialize()
    cmd.fragment("gly")
    r = psico.querying.iterate_to_list("not hydro", "name")
    assert r == ["N", "CA", "C", "O"]


def test_csp():
    cmd.reinitialize()
    cmd.fab("A/1/ ARKA B/1/ GERD", "m1")
    cmd.fab("EE", "m2")
    cmd.fab("DD", "m3")
    assert psico.querying.csp("m1") == -2
    assert psico.querying.csp("m2", "m3") == 4


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
