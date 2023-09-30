import psico.load_mtz_cctbx as m
import importlib.util
import pytest
from pathlib import Path
from pymol import cmd

DATA_PATH = Path(__file__).resolve().parent / 'data'


def test_load_mtz_cctbx():
    cmd.reinitialize()

    if importlib.util.find_spec("iotbx") is None:
        pytest.skip("cctbx not found")

    m.load_mtz_cctbx(DATA_PATH / "5e5z_phases.mtz")

    assert set(cmd.get_names()) == {
        '5e5z_phases_FC_ALL_LS_PHIC_ALL_LS',
        '5e5z_phases_DELFWT_PHDELWT',
        '5e5z_phases_FWT_PHWT',
        '5e5z_phases_FC_ALL_PHIC_ALL',
        '5e5z_phases_FC_PHIC',
    }
