import psico.moe
import psico.querying
from pathlib import Path
from pytest import approx
from pymol import cmd

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_load_moe():
    cmd.reinitialize()
    psico.moe.load_moe(DATA_PATH / "4ins-frag.moe", "m1")
    assert cmd.get_chains() == ["4INS.A_0", "4INS.B_1", "4INS.B_2", "4INS.B_3"]
    assert cmd.count_atoms() == 152
    assert cmd.count_atoms("solvent") == 8
    assert cmd.count_atoms("resn HOH") == 8
    assert cmd.count_atoms("inorganic") == 1
    assert cmd.count_atoms("elem Zn") == 1
    assert cmd.count_atoms("name Zn") == 1
    assert cmd.count_atoms("resn CYS and guide") == 4
    assert cmd.count_atoms("name SG and bound_to name SG") == 4
    assert cmd.count_atoms("chain 4INS.B_1 and polymer") == 61

    symmetry = cmd.get_symmetry()
    assert symmetry[:6] == approx([82.5, 82.5, 34, 90, 90, 120])
    assert symmetry[6] == "H3"

    assert cmd.get_names("selections") == ["set_SIC", "set_cys_sidechain"]
    resn_list = psico.querying.iterate_to_list("set_SIC", "resn")
    assert resn_list == ["SER"] * 6 + ["ILE"] * 8 + ["CYS"] * 6
    name_list = psico.querying.iterate_to_list("set_cys_sidechain", "name")
    assert name_list == ["CA", "CB", "SG"]
