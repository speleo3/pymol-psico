import psico.prolif
from pymol import cmd
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_prolif_interactions():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "2x19-frag-mse.pdb", "m1")
    cmd.pseudoatom("m2", elem="Mg", pos=(-27, -2, 6))
    psico.prolif.prolif_interactions("m2", "m1")
    names = cmd.get_names()
    assert "prolif.VdWContact" in names
