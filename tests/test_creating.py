import psico.creating
import psico.querying
import os
from pymol import cmd
from pytest import approx

DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')
FILENAME_WITH_MSE = os.path.join(DATA_DIR, '2x19-frag-mse.pdb')


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


# def join_states(name, selection='all', discrete=-1, zoom=0, quiet=1):
# def ramp_levels(name, levels, quiet=1):
# def pdb2pqr(name, selection='all', ff='amber', debump=1, opt=1, assignonly=0,
# def corina(name, selection, exe='corina', state=-1, preserve=0, quiet=1):
# def prepwizard(name, selection='all', options='', state=-1,
# def fiber(seq, num=4, name='', rna=0, single=0, repeats=0,
