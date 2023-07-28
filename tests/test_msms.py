import psico.msms
import psico.querying
import pytest
from pymol import cmd
from pytest import approx
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parent / 'data'


@pytest.mark.exe
def test_msms_surface():
    eps = 1e-3
    cmd.reinitialize()
    cmd.fragment('ala', 'm1')
    # default
    psico.msms.msms_surface(name='surf1')
    extent = cmd.get_extent('surf1')
    assert extent[0] == approx([-2.705, -3.208, -2.413], rel=eps)
    assert extent[1] == approx([3.530, 2.907, 2.676], rel=eps)
    # global solvent_radius
    cmd.set('solvent_radius', 3.5)
    psico.msms.msms_surface(name='surf2')
    extent = cmd.get_extent('surf2')
    assert extent[0] == approx([-2.705, -3.169, -2.436], rel=eps)
    assert extent[1] == approx([3.530, 2.907, 2.676], rel=eps)
    # object-level solvent_radius
    cmd.set('solvent_radius', 2.8, 'm1')
    psico.msms.msms_surface(name='surf3')
    extent = cmd.get_extent('surf3')
    assert extent[0] == approx([-2.705, -3.161, -2.427], rel=eps)
    assert extent[1] == approx([3.530, 2.907, 2.676], rel=eps)
    # modified atom radii
    cmd.alter('m1', 'vdw = 3.0')
    psico.msms.msms_surface(name='surf4')
    extent = cmd.get_extent('surf4')
    assert extent[0] == approx([-4.605, -5.162, -4.418], rel=eps)
    assert extent[1] == approx([5.030, 4.861, 4.681], rel=eps)

def test_save_xyzrn(tmp_path):
    cmd.reinitialize()
    cmd.fragment("gly")
    path = tmp_path / "f.xyzrn"
    colors = []
    psico.msms.save_xyzr(str(path), "not hydro", _colorsout=colors)
    assert path.read_text().splitlines() == [
        "-1.195 0.201 -0.206 1.80 1 N_GLY_22",
        "0.230 0.318 -0.502 1.80 1 CA_GLY_22",
        "1.059 -0.390 0.542 1.80 1 C_GLY_22",
        "0.545 -0.975 1.499 1.50 1 O_GLY_22",
    ]
    assert colors == [27, 26, 26, 28]


def test_atmtypenumbers():
    cmd.reinitialize()
    cmd.fab("AEV", "m1", hydro=0)
    psico.msms.atmtypenumbers(DATA_PATH / "atmtypenumbers", "m1")
    vdw = dict(psico.querying.iterate_to_list("m1", "(resn, name), vdw"))
    assert vdw["ALA", "CB"] == approx(2.0)
    assert vdw["GLU", "CB"] == approx(2.0)
    assert vdw["GLU", "OE1"] == approx(1.4)
    assert vdw["VAL", "CB"] == approx(2.0)
    psico.msms.atmtypenumbers(DATA_PATH / "atmtypenumbers", "m1", united=0)
    vdw = dict(psico.querying.iterate_to_list("m1", "(resn, name), vdw"))
    assert vdw["ALA", "CB"] == approx(1.74)
    assert vdw["GLU", "CB"] == approx(1.74)
    assert vdw["GLU", "OE1"] == approx(1.4)
    assert vdw["VAL", "CB"] == approx(1.74)
