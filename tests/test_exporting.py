import psico.exporting
from pymol import cmd


def test_save_traj(tmp_path):
    cmd.reinitialize()
    cmd.fragment('gly', 'm1')
    for state in range(1, 4):
        cmd.create('m2', 'm1', 1, state)

    for ext in ['dcd', 'crd']:
        filename = str(tmp_path / "out.") + ext
        # export
        psico.exporting.save_traj(filename, 'm2')
        # import (verify export)
        cmd.create('m3', 'm1')
        cmd.load_traj(filename, 'm3', 1)
        assert cmd.count_states('m3', 3)
        cmd.delete('m3')


# def get_pdb_sss(selection='(all)', state=-1, quiet=1):
# def get_pdb_seqres(selection='all', quiet=1):
# def save_pdb(filename, selection='(all)', state=-1, symm=1, ss=1, aniso=0, seqres=0, quiet=1):
# def save(filename, selection='(all)', state=-1, format='',
# def unittouu(string, dpi=90.0):
# def paper_png(filename, width=100, height=0, dpi=300, ray=1):
# def save_pdb_without_ter(filename, selection, *args, **kwargs):
