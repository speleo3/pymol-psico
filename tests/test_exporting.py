import psico.exporting
import pytest
from pytest import approx
from pymol import cmd

save_traj_params = [
    ('dcd', psico.exporting.save_traj),
    ('crd', psico.exporting.save_traj),
]

try:
    import mdtraj
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


# def get_pdb_sss(selection='(all)', state=-1, quiet=1):
# def get_pdb_seqres(selection='all', quiet=1):
# def save_pdb(filename, selection='(all)', state=-1, symm=1, ss=1, aniso=0, seqres=0, quiet=1):
# def save(filename, selection='(all)', state=-1, format='',
# def unittouu(string, dpi=90.0):
# def paper_png(filename, width=100, height=0, dpi=300, ray=1):
# def save_pdb_without_ter(filename, selection, *args, **kwargs):
