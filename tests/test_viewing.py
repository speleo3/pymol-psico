import psico.viewing
import psico.moe
from psico.querying import iterate_to_list
from pathlib import Path
from pymol import cmd

DATA_PATH = Path(__file__).resolve().parent / 'data'


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
