import psico.wizards
from pymol import cmd
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_names_sc():
    cmd.reinitialize()
    sc = psico.wizards.names_sc()
    assert sc.has_key("sspick")
    assert sc.has_key("message")


def test_Sspick():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1ubq.cif.gz", "m1")
    wiz = psico.wizards.Sspick()
    cmd.edit("guide & resi 28")
    wiz.do_pick(bondFlag=0)
    assert cmd.count_atoms("ss_01 & guide") == 12
    assert cmd.count_atoms("?pk1") == 0
