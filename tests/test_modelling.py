import psico.modelling
from pymol import cmd
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_mutate():
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')
    psico.modelling.mutate("resi 134", "K", inplace=1)
    assert psico.modelling.get_seq("all") == "LSKMPD"


def test_mutate_all():
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')
    psico.modelling.mutate_all("resi 134+137", "K", inplace=1)
    assert psico.modelling.get_seq("all") == "LSKMPK"


def test_add_missing_atoms():
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')
    cmd.remove("not backbone")
    assert cmd.count_atoms("resi 132") == 4
    psico.modelling.add_missing_atoms('resi 132+133', cycles=10)
    assert cmd.count_atoms("resi 132") == 8


def test_peptide_rebuild():
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb', "m1")
    psico.modelling.peptide_rebuild("m2", "m1", cycles=100)
    assert psico.modelling.get_seq("m2") == "LSMMPD"


def test_get_seq():
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')
    assert psico.modelling.get_seq("all") == "LSMMPD"
    cmd.remove("resi 134")
    assert psico.modelling.get_seq("all") == "LS/MPD"
    cmd.alter("resi 137", "resn = 'XXX'")
    assert psico.modelling.get_seq("all", unknown="#") == "LS/MP#"
