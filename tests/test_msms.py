import psico.msms
import pytest
from pymol import cmd
from pytest import approx


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
