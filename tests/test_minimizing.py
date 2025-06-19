import psico.minimizing
from pymol import cmd
import pytest
from pytest import approx


def test_get_fixed_indices():
    cmd.reinitialize()
    cmd.fragment("gly", "m1")
    cmd.flag("fix", "elem N+O")
    r = psico.minimizing.get_fixed_indices("m1", 1, _self=cmd)
    assert r == [0, 3]
    r = psico.minimizing.get_fixed_indices("m1 & elem C+O", 1, _self=cmd)
    assert r == [2]


def test_load_or_update():
    cmd.reinitialize()
    cmd.fragment("gly", "m1")
    molstr = cmd.get_str("mol", "m1")
    psico.minimizing.load_or_update(molstr, "m2", "m1", 1, _self=cmd)
    assert cmd.rms_cur("m2", "m1", matchmaker=-1) == approx(0.0)
    cmd.alter_state(1, "m1", "(x,y,z) = (x*2, y*2, z*2)")
    assert cmd.rms("m2", "m1", matchmaker=-1) > 1
    psico.minimizing.load_or_update(molstr, "", "m1", 1, _self=cmd)
    assert cmd.get_names() == ["m1", "m2"]
    assert cmd.rms("m2", "m1", matchmaker=-1) == approx(0.0)


@pytest.mark.openbabel
def test_minimize_ob():
    cmd.reinitialize()
    cmd.fab("ACD", "m1")
    cmd.copy("m2", "m1")
    psico.minimizing.minimize_ob("m2", nsteps=100, ff="UFF")
    r = cmd.rms("m2", "m1")
    assert r == approx(0.1739, abs=1e-2)
    psico.minimizing.minimize_ob("m2", nsteps=100, ff="MMFF94s")
    r = cmd.rms("m2", "m1")
    assert r == approx(0.2198, abs=1e-2)


@pytest.mark.rdkit
def test_minimize_rdkit():
    cmd.reinitialize()
    cmd.fab("ACD", "m1")
    cmd.copy("m2", "m1")
    psico.minimizing.minimize_rdkit("m2", nsteps=100, ff="MMFF94")
    r = cmd.rms("m2", "m1")
    assert r == approx(0.7118, abs=1e-2)
    psico.minimizing.minimize_rdkit("m2", nsteps=100, ff="UFF")
    r = cmd.rms("m2", "m1")
    abs_UFF = 0.03
    assert r == approx(0.9508, abs=abs_UFF)


@pytest.mark.openbabel
def test_clean_ob():
    cmd.reinitialize()
    cmd.fab("ACD", "m1")
    cmd.copy("m2", "m1")
    psico.minimizing.clean_ob("m2 & resn CYS+ASP", "m2")
    r = cmd.rms("m2 & resn ALA & not hydro", "m1 & resn ALA")
    assert r == approx(0.0)
    r = cmd.rms("m2", "m1")
    assert r == approx(0.1108, abs=1e-2)
