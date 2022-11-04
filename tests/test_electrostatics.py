import psico.electrostatics
from pymol import cmd
from pathlib import Path
import pytest

DATA_PATH = Path(__file__).resolve().parent / 'data'


@pytest.mark.exe
def test_apbs_surface():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "2x19-frag-mse.pdb", "m1")
    psico.electrostatics.apbs_surface(map_name="mapfoo", group_name="thegroup")
    names = cmd.get_names()
    assert "thegroup" in names
    assert "m1" in names
    assert "mapfoo" in names
    assert "ramp01" in names
    cmd.delete("thegroup")  # deletes the group and all its members
    assert cmd.get_names() == []
