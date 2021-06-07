import psico.modelling
from pymol import cmd
from pathlib import Path
from pytest import approx, mark

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


def test_update_align():
    cmd.reinitialize()
    cmd.fab("CDEGGGGGGKRC", "m1")
    cmd.fab("CDEAAKRC", "m2")
    psico.modelling.update_align("m1", "m2")
    assert (cmd.get_coords("/m1///1-3") == cmd.get_coords("/m2///1-3")).all()
    assert (cmd.get_coords("/m1///10-") == cmd.get_coords("/m2///6-8")).all()
    assert 6 < cmd.get_distance("/m1///6/N", "/m1///5/C")


@mark.parametrize("fix", ["fix", "protect"])
def test_sculpt_homolog__fix(fix):
    cmd.reinitialize()
    cmd.fab("CDEGGGGGGKRC", "m1")
    cmd.fab("CDEAAKRC", "m2")
    psico.modelling.sculpt_homolog("m1", "m2", cycles=10, fix=fix)
    assert (cmd.get_coords("/m1///1-3") == cmd.get_coords("/m2///1-3")).all()
    assert (cmd.get_coords("/m1///10-") == cmd.get_coords("/m2///6-8")).all()
    assert 6 > cmd.get_distance("/m1///6/N", "/m1///5/C")


def test_sculpt_homolog__restrain():
    cmd.reinitialize()
    cmd.fab("CDEGGGGGGKRC", "m1")
    cmd.fab("CDEAAKRC", "m2")
    psico.modelling.sculpt_homolog("m1", "m2", cycles=100, fix="restrain")
    c1 = cmd.get_coords("/m1///10-")
    c2 = cmd.get_coords("/m2///6-8")
    assert not (c1 == c2).all()
    assert c1.flatten().tolist() == approx(c2.flatten().tolist(), abs=1.5)
