import psico.exporting
import pytest
from pytest import approx
from pymol import cmd
from pathlib import Path
from typing import List, Tuple, Callable

DATA_PATH = Path(__file__).resolve().parent / 'data'

save_traj_params: List[Tuple[str, Callable[[str, str], object]]] = [
    ('dcd', psico.exporting.save_traj),
    ('crd', psico.exporting.save_traj),
]

try:
    import mdtraj  # noqa: F401
except ImportError:
    print("Skipping mdtraj tests")
else:
    save_traj_params += [
        ('dcd', psico.exporting.save_mdtraj),
        ('xtc', psico.exporting.save_mdtraj),
    ]


@pytest.mark.parametrize("ext,save_func", save_traj_params)
def test_save_traj(ext, save_func, tmp_path):
    cmd.reinitialize()
    cmd.fragment('gly', 'm1')
    for state in range(1, 4):
        cmd.create('m2', 'm1', 1, state)
        cmd.rotate("y", 30 * state, "m2", state)

    filename = str(tmp_path / "out.") + ext
    # export
    save_func(filename, 'm2')
    # import (verify export)
    cmd.create('m3', 'm1')
    cmd.load_traj(filename, 'm3', 1)
    assert cmd.count_states('m3', 3)
    coords3 = cmd.get_coords('m3', 0)
    coords2 = cmd.get_coords('m2', 0)
    assert coords3.shape == coords2.shape
    assert coords3 == approx(coords2, abs=1e-2)


def test_save_pdb(tmp_path):
    cmd.reinitialize()
    cmd.load(DATA_PATH / '1ubq.cif.gz')
    filename = tmp_path / "out.pdb"
    psico.exporting.save_pdb(filename, symm=1, ss=1, seqres=1)
    content = filename.read_text()
    assert "\nSEQRES   2 A   76  THR LEU GLU VAL" in content
    assert "\nSHEET    1   1 1 MET A   1  THR A   7  0" in content
    assert "\nCRYST1   50.840   42.770   28.950  90.00  90.00  90.00 P 21 21 21" in content


def test_save_pdb_without_ter(tmp_path):
    cmd.reinitialize()
    cmd.load(DATA_PATH / '2x19-frag-mse.pdb')
    cmd.remove("resi 135")
    cmd.alter("resi 136-137", "chain='C'")
    assert cmd.get_pdbstr().count("TER") == 2
    filename = tmp_path / "pdbwithoutter.pdb"
    psico.exporting.save_pdb_without_ter(filename, "all")
    with open(filename) as handle:
        content = handle.read()
    assert content.count("TER") == 0
    assert content.count("ATOM") == 26


@pytest.mark.xfail
def test_get_grostr():
    cmd.reinitialize()
    filename = DATA_PATH / 'ala-cys-asp.gro'
    orig = filename.read_text()
    cmd.load(filename, "m1")
    grostr = psico.exporting.get_grostr()
    # ignore first line (title)
    assert orig.splitlines()[1:] == grostr.splitlines()[1:]


def test_get_grostr__dims_from_extent():
    cmd.reinitialize()
    cmd.fragment("gly")
    grostr = psico.exporting.get_grostr()
    assert grostr.splitlines()[-1] == "   0.26170   0.23123   0.29780"


def test_unittouu():
    unittouu = psico.exporting.unittouu
    assert unittouu("1in") == approx(90)
    assert unittouu("2.54cm") == approx(90)
    assert unittouu("25.4mm", dpi=150) == approx(150)
    with pytest.raises(ValueError) as excinfo:
        unittouu(None)
    assert 'cannot parse value' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        unittouu("#")
    assert 'cannot parse value' in str(excinfo.value)
    with pytest.raises(ValueError) as excinfo:
        unittouu("123xx")
    assert 'unknown unit' in str(excinfo.value)


# def get_pdb_sss(selection='(all)', state=-1, quiet=1):
# def get_pdb_seqres(selection='all', quiet=1):
# def save(filename, selection='(all)', state=-1, format='',
# def paper_png(filename, width=100, height=0, dpi=300, ray=1):
