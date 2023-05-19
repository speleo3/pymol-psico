import psico.viewing
import psico.moe
import pymol
import pytest
from pytest import approx
from psico.querying import iterate_to_list
from pathlib import Path
from pymol import cmd

DATA_PATH = Path(__file__).resolve().parent / 'data'
FILENAME_WITH_MSE = DATA_PATH / '2x19-frag-mse.pdb'


def test_nice():
    cmd.reinitialize()
    cmd.load(FILENAME_WITH_MSE)
    psico.viewing.nice()
    assert cmd.count_atoms("rep cartoon") > 0
    assert cmd.count_atoms("rep ribbon") == 0
    psico.viewing.nice(simple=10)
    assert cmd.count_atoms("rep cartoon") == 0
    assert cmd.count_atoms("rep ribbon") > 0


def test_get_color_family():
    blues = psico.viewing.get_color_family("blues")
    assert "blue" in blues
    assert "marine" in blues
    assert "red" not in blues
    with pytest.raises(pymol.CmdException):
        psico.viewing.get_color_family("nosuchcolors")


def tst_cbm():
    cmd.reinitialize()
    cmd.fragment("gly", "m1")
    cmd.fragment("gly", "m2")
    cmd.fragment("gly", "m3")
    psico.viewing.cbm("*")
    assert set(iterate_to_list("*", "color")) == {2, 3, 4}
    psico.viewing.cbm("*", "reds")
    assert set(iterate_to_list("*", "color")) == {4, 32, 5268}


def tst_cbs():
    cmd.reinitialize()
    cmd.fab("ACDEFG")
    cmd.alter("resi 1-3", "segi = 'A'")
    cmd.alter("resi 3-6", "segi = 'B'")
    psico.viewing.cbs("*")
    assert set(iterate_to_list("*", "color")) == {2, 3}
    psico.viewing.cbs("*", "reds")
    assert set(iterate_to_list("*", "color")) == {4, 32}


def test_spectrumany():
    cmd.reinitialize()
    cmd.fab("ACDEFG")
    psico.viewing.spectrumany("resi", "blue black")
    colors = iterate_to_list("(first resi 1)(last resi 6)", "color")
    assert colors == approx([0x400000FF, 0x40000000], abs=60)  # very fuzzy
    colors = iterate_to_list("(first resi 2)(last resi 5)", "color")
    assert colors != approx([0x400000FF, 0x40000000], abs=60)  # very fuzzy
    psico.viewing.spectrumany("resi", "blue black", minimum=2, maximum=5)
    colors = iterate_to_list("(first resi 2)(last resi 5)", "color")
    assert colors == approx([0x400000FF, 0x40000000], abs=60)  # very fuzzy


def test_spectrum_states():
    cmd.reinitialize()
    cmd.pseudoatom("m1")
    cmd.create("m1", "m1", 1, 2)
    cmd.create("m1", "m1", 1, 3)
    psico.viewing.spectrum_states(representations="spheres",
                                  color_list="blue red yellow")
    assert cmd.get_setting_int("sphere_color", "m1", 1) == 0x400000FF
    assert cmd.get_setting_int("sphere_color", "m1", 2) == 0x40FF0000
    assert cmd.get_setting_int("sphere_color", "m1", 3) == 0x40FFFF00


def test_scene_preserve():
    cmd.reinitialize()
    cmd.pseudoatom("m1")
    cmd.show_as("spheres")
    with psico.viewing.scene_preserve():
        cmd.hide("everything")
        assert cmd.count_atoms("rep spheres") == 0
    assert cmd.count_atoms("rep spheres") == 1


def test_goodsell_lighting():
    cmd.reinitialize()
    psico.viewing.goodsell_lighting()
    assert cmd.get_setting_float("ambient") == approx(0.6)
    assert cmd.get_setting_int("ray_shadow") == 0
    assert cmd.get_setting_float("specular_intensity") == approx(0)


def test_grid_labels():
    cmd.reinitialize()
    cmd.fragment("ala")
    cmd.fragment("gly")
    cmd.remove("hydro")
    psico.viewing.grid_labels()
    labels = iterate_to_list("all", "label")
    assert labels == ["ala", "", "", "", "", "gly", "", "", ""]
    cmd.label()
    cmd.set_title("ala", 1, "foo")
    cmd.set_title("gly", 1, "bar")
    psico.viewing.grid_labels("title", subsele="name CA")
    labels = iterate_to_list("all", "label")
    assert labels == ["", "foo", "", "", "", "", "bar", "", ""]
    cmd.label()
    cmd.group("g1", "ala gly")
    psico.viewing.grid_labels()
    labels = iterate_to_list("all", "label")
    assert labels == ["ala", "", "", "", "", "", "", "", ""]
    cmd.label()
    psico.viewing.grid_labels("'abc'", "gly")
    labels = iterate_to_list("all", "label")
    assert labels == ["", "", "", "", "", "abc", "", "", ""]


def test_show_ptm():
    cmd.reinitialize()
    psico.moe.load_moe(DATA_PATH / "4ins-frag.moe", "m1")
    cmd.hide()
    psico.viewing.show_ptm()
    assert cmd.count_atoms("rep spheres") == 4
    assert cmd.count_atoms("rep sticks") == 24
