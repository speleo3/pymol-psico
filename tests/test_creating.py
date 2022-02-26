import psico.creating
import psico.querying
from pathlib import Path
from pymol import cmd
from pytest import approx

DATA_PATH = Path(__file__).resolve().parent / 'data'
FILENAME_WITH_MSE = DATA_PATH / '2x19-frag-mse.pdb'


def test_sidechaincenters():
    cmd.reinitialize()
    cmd.load(FILENAME_WITH_MSE)
    cmd.alter("all", "b = 5")
    cmd.alter("resn SER", "b = 20")
    psico.creating.sidechaincenters(name='SCA')
    psico.creating.sidechaincenters("scc_centroid", method="centroid")
    assert cmd.count_atoms("scc") == 5
    assert cmd.count_atoms("scc & resn MSE") == 2
    assert cmd.count_atoms("scc & resn ASP") == 0
    assert cmd.count_atoms("name SCA") == 5
    assert cmd.count_atoms("scc_centroid") == 6
    assert cmd.count_atoms("scc_centroid & resn MSE") == 2
    assert cmd.count_atoms("scc_centroid & resn ASP") == 1
    assert cmd.count_atoms("name PS1") == 6
    r = psico.querying.iterate_to_list("scc", "b")
    assert r == approx([5, 20, 5, 5, 5])


def test_join_states():
    cmd.reinitialize()
    cmd.load(DATA_PATH / "1nmr-frag-nohydro.pdb", "m1", multiplex=1)
    assert len(cmd.get_names()) == 4
    assert cmd.count_states() == 1

    psico.creating.join_states("m2", "m1_*", discrete=0)
    assert cmd.count_states("m2") == 4
    assert cmd.count_discrete("m2") == 0
    assert cmd.count_atoms("m2") == 241

    psico.creating.join_states("m3", "m1_*", discrete=1)
    assert cmd.count_states("m3") == 4
    assert cmd.count_discrete("m3") == 1
    assert cmd.count_atoms("m3") == 241 * 4


# def ramp_levels(name, levels, quiet=1):
# def pdb2pqr(name, selection='all', ff='amber', debump=1, opt=1, assignonly=0,
# def corina(name, selection, exe='corina', state=-1, preserve=0, quiet=1):
# def prepwizard(name, selection='all', options='', state=-1,
# def fiber(seq, num=4, name='', rna=0, single=0, repeats=0,
