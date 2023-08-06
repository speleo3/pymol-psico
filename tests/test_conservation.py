import psico.conservation
import psico.querying
from pathlib import Path
from pymol import cmd
import pytest

DATA_PATH = Path(__file__).resolve().parent / 'data'


@pytest.mark.web
def test_consurfdb():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1ubq.cif.gz", "m1")
    cmd.color("white")
    psico.conservation.consurfdb("6Q00", "A", "m1")
    colors = psico.querying.iterate_to_list("m1", "color")
    assert len(set(colors)) > 10


def test_load_consurf():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1ubq.cif.gz", "m1")
    cmd.color("white")
    psico.conservation.load_consurf(
        str(DATA_PATH / "consurf_summary-5NVGA.txt"), "m1")
    colors = psico.querying.iterate_to_list("m1", "color")
    assert len(set(colors)) > 10
