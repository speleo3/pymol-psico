import psico.morphing
from pathlib import Path

DATA_PATH = Path(__file__).resolve().parent / "data"


def test_morpheasy_linear():
    from pymol2 import PyMOL
    with PyMOL() as p1:
        p1.cmd.load(DATA_PATH / "1nmr-frag-nohydro.pdb", "m1")
        p1.cmd.create("m2", "m1", 2, 1)
        psico.morphing.morpheasy_linear("m1", "m2",
                                        source_state=1,
                                        target_state=1,
                                        name="m3", _self=p1.cmd)
        assert p1.cmd.count_states("m3") == 30
