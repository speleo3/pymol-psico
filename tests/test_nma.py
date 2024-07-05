import psico.nma
from pymol import cmd


def test_normalmodes_prody():
    cmd.reinitialize()
    cmd.fab("ACDEFGH", "m1", ss=1)
    psico.nma.normalmodes_prody("m1", prefix="pro", states=4)
    assert cmd.get_names() == ["m1", "pro7", "pro8", "pro9", "pro10"]
    assert cmd.count_states() == 4
    assert cmd.count_atoms("pro7") == 7
