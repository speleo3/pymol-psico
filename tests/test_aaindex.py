import psico.aaindex
import psico.querying
from pymol import cmd
from pathlib import Path
from pytest import approx

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_hydropathy2b(tmp_path):
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')

    # use cached sample file, skip download
    cmd.set("fetch_path", str(DATA_PATH / "aaindex"))

    psico.aaindex.hydropathy2b()

    b_list = psico.querying.iterate_to_list("guide", "b")
    assert b_list == approx([3.8, -0.8, 1.9, 1.9, -1.6, -3.5])


def test_aaindex2b(tmp_path):
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')

    # use cached sample file, skip download
    cmd.set("fetch_path", str(DATA_PATH / "aaindex"))

    psico.aaindex.aaindex2b("LAWE840101", "name CA", "q")

    b_list = psico.querying.iterate_to_list("guide", "q")
    assert b_list == approx([1.02, 0.05, 0.81, 0.81, 2.03, -0.75])
