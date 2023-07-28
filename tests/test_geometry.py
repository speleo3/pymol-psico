import psico.geometry
from pymol import cmd
import pytest
from pytest import approx


@pytest.mark.exe
def test_delaunay():
    cmd.reinitialize()
    cmd.fab("GGGGG", "m1", ss=1)
    psico.geometry.delaunay("guide & m1", name="g1")
    assert cmd.count_atoms("g1") == 5
    assert cmd.get_extent("g1") == [
        approx([0.230, -4.570, -1.174], abs=1e-3),
        approx([5.335, 0.318, 2.799], abs=1e-3),
    ]
    psico.geometry.delaunay("guide & m1", name="c1", as_cgo=1)
    assert cmd.get_type("c1") == "object:cgo"
    assert cmd.get_extent("c1") == [
        approx([0.180, -4.620, -1.224], abs=1e-3),
        approx([5.385, 0.368, 2.849], abs=1e-3),
    ]
